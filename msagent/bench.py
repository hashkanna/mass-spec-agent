"""Retrieval benchmark: natural-language questions with explicitly defined gold studies.

Each question in benchmarks/queries.yaml carries a ``gold_sql`` rule written by hand against the raw
SDRF values (not through msagent's vocabulary), so the rule is inspectable and editable. Systems:

* ``bm25``      keyword search over the same per-study metadata documents (stand-in for keyword search)
* ``semantic``  embedding similarity only ("LLM/vector search in front of the archive")
* ``hybrid``    BM25 + embeddings, reciprocal-rank fusion
* ``msagent``   criteria parsing + facet graph + sample-level verification (this project)
* ``pride``     PRIDE's own search API, only with --pride-live and network access to www.ebi.ac.uk

recall@k is hits in the top k divided by min(k, |gold|), so 1.0 is reachable when the gold set is
larger than k (a normalised variant). standard_recall@10 is hits in the top 10 divided by |gold|. precision@10 is the share of the (up to) 10 returned studies that are gold.
"""

from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

from .cohort import find_cohorts
from .store import Store

# Short views over the raw tables so gold rules stay readable.
GOLD_PRELUDE = """
WITH sa AS (SELECT accession, source_name, attribute AS a, lower(value) AS v FROM sample_attributes),
runs AS (SELECT accession, source_name, lower(coalesce(instrument, '')) AS instrument, lower(coalesce(label, '')) AS label,
                lower(coalesce(acquisition, '')) AS acq, file_format FROM assays),
mods AS (SELECT accession, lower(name) AS name, upper(coalesce(targets, '')) AS targets, unimod FROM modifications),
proj AS (SELECT * FROM projects)
"""

SYSTEMS = [
    ("bm25", "Keyword (BM25)"),
    ("semantic", "Semantic (embeddings)"),
    ("hybrid", "Hybrid keyword+semantic"),
    ("msagent", "Mass-Spec Agent"),
]
KS = (10, 50)


def gold_set(store: Store, sql: str) -> set[str]:
    return {r[0] for r in store.sql(GOLD_PRELUDE + sql)}


def score(ranking: list[str], gold: set[str]) -> dict:
    out = {}
    for k in KS:
        hits = sum(1 for a in ranking[:k] if a in gold)
        out[f"recall@{k}"] = round(hits / min(k, len(gold)), 4) if gold else 0.0
    top = ranking[:10]
    out["standard_recall@10"] = round(sum(1 for a in ranking[:10] if a in gold) / len(gold), 4) if gold else 0.0
    out["precision@10"] = round(sum(1 for a in top if a in gold) / len(top), 4) if top else 0.0
    first = next((i for i, a in enumerate(ranking) if a in gold), None)
    out["mrr"] = round(1 / (first + 1), 4) if first is not None else 0.0
    return out


def rank_system(store: Store, system: str, question: str) -> list[str]:
    if system == "bm25":
        return [a for a, _ in store.text.keyword(question, 50)]
    if system == "semantic":
        return [a for a, _ in store.text.semantic(question, 50)]
    if system == "hybrid":
        return [a for a, _ in store.text.hybrid(question, 50)]
    if system == "msagent":
        res = find_cohorts(store, question, limit=50, verify_n=80)
        return [r["accession"] for r in res["results"]]
    if system == "pride":
        from .pride_api import search_projects

        return search_projects(question, page_size=50)
    raise ValueError(system)


def run(queries_path: str, output_path: str, pride_live: bool = False) -> dict:
    store = Store()
    spec = yaml.safe_load(Path(queries_path).read_text())
    systems = list(SYSTEMS) + ([("pride", "PRIDE search API (live)")] if pride_live else [])
    rows = []
    for q in spec["queries"]:
        gold = gold_set(store, q["gold_sql"])
        if q.get("must_include"):
            gold |= set(q["must_include"])
        entry = {"id": q["id"], "question": q["question"], "gold_rule": q.get("gold_rule", ""), "n_gold": len(gold), "scores": {}, "top": {}, "seconds": {}}
        for sid, _ in systems:
            t0 = time.time()
            ranking = rank_system(store, sid, q["question"])
            entry["seconds"][sid] = round(time.time() - t0, 2)
            entry["scores"][sid] = score(ranking, gold)
            entry["top"][sid] = ranking[:10]
        rows.append(entry)
        print(f"{q['id']:<28} gold={len(gold):<4} " + "  ".join(
            f"{sid}:{entry['scores'][sid]['recall@10']:.2f}/{entry['scores'][sid]['recall@50']:.2f}" for sid, _ in systems), flush=True)
    scored = [r for r in rows if r["n_gold"]]
    aggregate = {
        sid: {m: round(sum(r["scores"][sid][m] for r in scored) / max(1, len(scored)), 4) for m in ("recall@10", "standard_recall@10", "recall@50", "precision@10", "mrr")}
        for sid, _ in systems
    }
    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": spec.get("note", ""),
        "index": {"sdrf_commit": store.meta.get("sdrf_commit"), "embedder": store.text.meta.get("embedder")},
        "systems": [{"id": s, "label": label} for s, label in systems],
        "aggregate": aggregate,
        "queries": rows,
    }
    Path(output_path).write_text(json.dumps(result, indent=1))
    print("\nmean over", len(scored), "questions:")
    for sid, label in systems:
        a = aggregate[sid]
        print(f"  {label:<28} recall@10 {a['recall@10']:.2f}  std-recall@10 {a['standard_recall@10']:.2f}  recall@50 {a['recall@50']:.2f}  precision@10 {a['precision@10']:.2f}  MRR {a['mrr']:.2f}")
    return result

"""Export a selected cohort as a sample/file manifest with provenance."""

from __future__ import annotations

import csv
import io
import json
from datetime import datetime, timezone

from .criteria import Criteria, parse_query
from .project import _url_for
from .store import Store
from .verify import verify

KEY_ATTRS = ["organism", "organism part", "disease", "cell type", "cell line", "sex", "age", "treatment", "compound", "genotype"]
COLUMNS = ["accession", "group", "group_evidence", "individual", "sample", "assay_name", "data_file", "file_format",
           "fraction", "label", "technical_replicate", "instrument", "url", "url_source"] + [f"char:{a}" for a in KEY_ATTRS] + ["sdrf_file"]


def build_manifest(store: Store, accessions: list[str], criteria: Criteria | str | None = None) -> dict:
    c = criteria if isinstance(criteria, Criteria) else parse_query(criteria or "")
    data = store.load_projects(accessions)
    rows, studies = [], []
    for acc in accessions:
        pd = data.get(acc)
        if pd is None:
            studies.append({"accession": acc, "error": "not indexed"})
            continue
        v = verify(pd, c)
        chosen = {x["sample"]: x for x in v["cohort"]}
        rows_ok = set(v["cohort_rows"])
        for a in pd.assays:
            s = chosen.get(a["source_name"])
            if s is None or a["row_id"] not in rows_ok:
                continue
            url, src = _url_for(pd, a)
            attrs = pd.samples.get(a["source_name"], {})
            row = {
                "accession": acc, "group": s["group"] or "", "group_evidence": s["group_evidence"] or "",
                "individual": s["individual"] or "", "sample": a["source_name"], "assay_name": a["assay_name"],
                "data_file": a["data_file"], "file_format": a["file_format"], "fraction": a["fraction"], "label": a["label"],
                "technical_replicate": a["technical_replicate"], "instrument": a["instrument"], "url": url or "",
                "url_source": src, "sdrf_file": a["sdrf_file"],
            }
            for k in KEY_ATTRS:
                row[f"char:{k}"] = "; ".join(attrs.get(k, []))
            rows.append(row)
        studies.append({
            "accession": acc, "status": v["status"], "design": v["design"], "checks": v["checks"],
            "counts": v["counts"], "links": store.links(acc, pd.info["sdrf_files"]),
        })
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "query": c.text,
        "criteria": c.to_dict(),
        "provenance": {
            "sdrf_source": store.meta.get("sdrf_source"),
            "sdrf_commit": store.meta.get("sdrf_commit"),
            "index_built_at": store.meta.get("built_at"),
        },
        "studies": studies,
        "n_rows": len(rows),
        "n_files": len({(r["accession"], r["data_file"]) for r in rows}),
        "rows": rows,
    }


def to_tsv(manifest: dict) -> str:
    buf = io.StringIO()
    buf.write(f"# msagent cohort manifest generated {manifest['generated_at']}\n")
    buf.write(f"# query: {manifest['query']}\n")
    buf.write(f"# SDRF source: {manifest['provenance']['sdrf_source']} @ {manifest['provenance']['sdrf_commit']}\n")
    for s in manifest["studies"]:
        if "error" in s:
            continue
        buf.write(f"# {s['accession']}: {s['status']} | {s['design']}\n")
        for ck in s["checks"]:
            buf.write(f"#   [{ck['status']}] {ck['label']}: {ck['detail']}\n")
    w = csv.DictWriter(buf, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for r in manifest["rows"]:
        w.writerow(r)
    return buf.getvalue()


def to_json(manifest: dict) -> str:
    return json.dumps(manifest, indent=1, default=str)

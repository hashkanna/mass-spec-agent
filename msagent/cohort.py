"""Find studies that can supply the requested cohort.

Pipeline:
  1. interpret the question into explicit criteria (``criteria.parse_query``),
  2. retrieve candidates from the facet graph (structured) and the text indexes (keyword + semantic),
  3. verify each shortlisted study criterion by criterion at sample level (``verify.verify``),
  4. rank: qualifying first, then studies needing review, then near misses with their reasons.
"""

from __future__ import annotations

import time
from collections import defaultdict

from . import vocab
from .criteria import CELL_TYPES, Criteria, parse_query
from .store import Store
from .verify import verify

STATUS_RANK = {"qualifies": 0, "needs review": 1, "excluded": 2}


def _values_matching(store: Store, facets: tuple[str, ...], pred) -> set[str]:
    out: set[str] = set()
    for f in facets:
        for value, accs in store.facets.get(f, {}).items():
            if pred(value):
                out.update(accs)
    return out


def _annotated(store: Store, facets: tuple[str, ...]) -> set[str]:
    out: set[str] = set()
    for f in facets:
        out |= store.annotated.get(f, set())
    return out


def facet_plan(store: Store, c: Criteria) -> list[tuple[str, set[str], set[str], float]]:
    """Project-level checks: (criterion id, passing accessions, annotated accessions, weight)."""
    plan = []
    if c.organisms:
        orgs = set(c.organisms)
        plan.append(("organism", _values_matching(store, ("organism",), lambda v: v in orgs), _annotated(store, ("organism",)), 1.0))
    if c.tissues:
        pats = [vocab.TISSUES.get(t) or vocab.TermPattern(t, [vocab.norm(t)]) for t in c.tissues]
        plan.append(("tissue", _values_matching(store, ("organism_part",), lambda v: any(p.matches(v) for p in pats)), _annotated(store, ("organism_part",)), 1.0))
    if c.diseases:
        pats = [vocab.DISEASES.get(d) or vocab.TermPattern(d, [vocab.norm(d)]) for d in c.diseases]
        plan.append(("disease", _values_matching(store, ("disease", "phenotype"), lambda v: any(p.matches(v) for p in pats)), _annotated(store, ("disease",)), 1.0))
    if c.comparison == "tumour_vs_normal" or c.comparison == "disease_vs_control":
        normal = _values_matching(store, ("disease", "phenotype", "sample_type", "organism_part"), vocab.is_normal_value)
        plan.append(("comparison", normal, _annotated(store, ("disease",)), 0.7))
    elif c.comparison == "treated_vs_control":
        ctrl = _values_matching(store, ("treatment", "compound"), lambda v: bool(vocab.CONTROL_RE.search(v)) or v in ("none", "no treatment"))
        treated = _values_matching(store, ("treatment", "compound"), lambda v: not vocab.CONTROL_RE.search(v))
        plan.append(("comparison", ctrl & treated, _annotated(store, ("treatment", "compound")), 0.7))
    elif c.comparison == "mutant_vs_wildtype":
        wt = _values_matching(store, ("genotype", "phenotype"), lambda v: bool(vocab.WILDTYPE_RE.search(v)))
        mut = _values_matching(store, ("genotype", "phenotype"), lambda v: not vocab.WILDTYPE_RE.search(v))
        plan.append(("comparison", wt & mut, _annotated(store, ("genotype",)), 0.7))
    if c.cell_types:
        pats = [CELL_TYPES.get(t) or vocab.TermPattern(t, [t]) for t in c.cell_types]
        plan.append(("cell_type", _values_matching(store, ("cell_type",), lambda v: any(p.matches(v) for p in pats)), _annotated(store, ("cell_type",)), 1.0))
    if c.cell_lines:
        pats = [vocab.TermPattern(x, [x]) for x in c.cell_lines]
        plan.append(("cell_line", _values_matching(store, ("cell_line", "cell_type"), lambda v: any(p.matches(v) for p in pats)), _annotated(store, ("cell_line",)), 1.0))
    if c.only_cell_lines:
        cl = _annotated(store, ("cell_line",)) | _values_matching(store, ("material_type", "organism_part"), lambda v: "cell" in v)
        plan.append(("only_cell_lines", cl, set(store.projects), 0.5))
    if c.instruments:
        subs = [x for i in c.instruments for x in vocab.INSTRUMENT_QUERY.get(i, [vocab.norm(i)])]
        plan.append(("instrument", _values_matching(store, ("instrument",), lambda v: any(s in vocab.norm(v) for s in subs)), _annotated(store, ("instrument",)), 1.0))
    if c.labels:
        labs = set(c.labels)
        ok = _values_matching(store, ("label_type",), lambda v: v in labs)
        if c.tmt_plex:
            ok &= _values_matching(store, ("tmt_plex",), lambda v: v == c.tmt_plex)
        plan.append(("label", ok, _annotated(store, ("label_type",)), 1.0))
    if c.acquisition:
        acq = set(c.acquisition)
        plan.append(("acquisition", _values_matching(store, ("acquisition",), lambda v: v in acq), _annotated(store, ("acquisition",)), 1.0))
    if c.ptms:
        want = (lambda v: True) if "any" in c.ptms else (lambda v: v in c.ptms)
        strong = _values_matching(store, ("ptm_enrichment",), want)
        weak = _values_matching(store, ("ptm",), want)
        plan.append(("ptm", strong, set(store.projects), 1.0))
        plan.append(("ptm_searched", weak, set(store.projects), 0.4))
    if c.raw_files:
        fmts = vocab.RAW_FORMATS | vocab.OPEN_FORMATS
        plan.append(("raw_files", _values_matching(store, ("file_format",), lambda v: v in fmts), _annotated(store, ("file_format",)), 0.5))
    if c.published_after:
        years = {a: (store.pride_years.get(a) or r[3]) for a, r in store.projects.items()}
        plan.append(("published_after", {a for a, y in years.items() if y and y >= c.published_after}, {a for a, y in years.items() if y}, 0.6))
    if c.min_individuals:
        rows = store.sql("SELECT accession, n_individuals, n_samples FROM projects")
        ok = {a for a, ni, ns in rows if ni is not None and ni >= c.min_individuals}
        known = {a for a, ni, ns in rows if ni is not None or ns < c.min_individuals}
        plan.append(("min_individuals", ok, known, 0.6))
    return plan


def find_cohorts(
    store: Store,
    query: str | Criteria,
    limit: int = 20,
    verify_n: int = 60,
    include_excluded: bool = True,
) -> dict:
    t0 = time.time()
    c = parse_query(query, [v for v, _ in store.facet_values("organism")]) if isinstance(query, str) else query
    text = c.text or " ".join(c.keywords)

    # ---- stage 1: project-level facet scoring + text relevance ---------------------------
    plan = facet_plan(store, c)
    score: dict[str, float] = defaultdict(float)
    max_score = sum(w for _, _, _, w in plan)
    for cid, ok, annotated, w in plan:
        for a in ok:
            score[a] += w
        for a in set(store.projects) - annotated:
            score[a] += 0.35 * w  # unknown at project level: keep in play, below confirmed
    text_rank: dict[str, float] = {}
    if store.text is not None and text.strip():
        for acc, s in store.text.hybrid(text, k=300):
            text_rank[acc] = s
    candidates = set(store.projects) if plan else set(text_rank)
    if plan:
        threshold = max_score - 1.0 - 1e-9 if len(plan) > 1 else max_score * 0.99
        candidates = {a for a in candidates if score.get(a, 0) >= threshold} | {a for a in text_rank if score.get(a, 0) >= max_score * 0.5}

    def stage1_key(a):
        return (-score.get(a, 0), -text_rank.get(a, 0), -store.projects[a][1])

    shortlist = sorted(candidates, key=stage1_key)[:verify_n]

    # ---- stage 2: sample-level verification -----------------------------------------------
    data = store.load_projects(shortlist)
    results = []
    for acc in shortlist:
        pd = data.get(acc)
        if pd is None:
            continue
        v = verify(pd, c)
        cohort = v.pop("cohort")
        results.append(
            {
                "accession": acc,
                "repository": pd.info["repository"],
                "title": (pd.pride or {}).get("title"),
                **v,
                "cohort_preview": cohort[:12],
                "text_score": round(text_rank.get(acc, 0.0), 5),
                "facet_score": round(score.get(acc, 0.0), 2),
                "links": store.links(acc, pd.info["sdrf_files"]),
                "sdrf_files": pd.info["sdrf_files"],
            }
        )

    def size(r):
        cnt = r["counts"]
        return cnt["individuals"] or cnt["cohort_samples"]

    def caveats(r):
        return sum(1 for ck in r["checks"] if ck.get("hint"))

    results.sort(key=lambda r: (STATUS_RANK[r["status"]], r["n_excluded"], r["n_unknown"], -r["n_confirmed"], caveats(r), -r["text_score"], -size(r)))
    kept = [r for r in results if r["status"] != "excluded"]
    near = [r for r in results if r["status"] == "excluded"]
    out = kept[:limit]
    if include_excluded and len(out) < limit:
        out += near[: max(3, limit - len(out))]
    return {
        "query": c.text,
        "criteria": c.to_dict(),
        "criteria_list": c.describe(),
        "summary": {
            "projects_indexed": len(store.projects),
            "facet_candidates": len(candidates),
            "verified": len(results),
            "qualifies": sum(1 for r in results if r["status"] == "qualifies"),
            "needs_review": sum(1 for r in results if r["status"] == "needs review"),
            "excluded": len(near),
            "seconds": round(time.time() - t0, 2),
        },
        "results": out,
    }


def graph_query(store: Store, limit: int = 50, min_samples: int | None = None, **filters) -> dict:
    """Deterministic facet filter over the metadata graph (all given filters must hold at project level).

    Filters take the same names/values as ``Criteria`` (organisms, tissues, diseases, instruments, labels,
    acquisition, ptms, cell_types, cell_lines, published_after, min_individuals).
    """
    c = Criteria.from_dict({k: v for k, v in filters.items() if v not in (None, [], "")})
    plan = [p for p in facet_plan(store, c) if p[0] not in ("ptm_searched", "comparison")]
    if not plan and not min_samples:
        return {"error": "give at least one filter", "filters": c.to_dict()}
    accs = set(store.projects)
    for _, ok, _, _ in plan:
        accs &= ok
    if min_samples:
        accs = {a for a in accs if store.projects[a][1] >= min_samples}
    rows = sorted(accs, key=lambda a: -store.projects[a][1])
    out = []
    for a in rows[:limit]:
        _, n_samples, n_runs, year = store.projects[a]
        out.append({"accession": a, "samples": n_samples, "runs": n_runs, "published": store.pride_years.get(a) or year})
    return {"filters": {k: v for k, v in c.to_dict().items() if v and k != "exclude_flagged"}, "n_matches": len(rows), "results": out}

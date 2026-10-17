"""Inspect one study: experimental design, and the sample -> run -> file map with filters.

Answers questions like "Which files in PXD012345 are the wild-type controls?".
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict

from . import vocab
from .store import Store
from .verify import classify_samples

# words in a filter that refer to an inferred role rather than a literal value
ROLE_WORDS = {
    "tumour": ("disease_role", "tumour"), "tumor": ("disease_role", "tumour"), "tumours": ("disease_role", "tumour"),
    "tumors": ("disease_role", "tumour"), "cancer": ("disease_role", "tumour"),
    "normal": ("disease_role", "normal"), "healthy": ("disease_role", "normal"), "non-tumour": ("disease_role", "normal"),
    "non-tumor": ("disease_role", "normal"), "adjacent": ("disease_role", "normal"),
    "treated": ("treatment_role", "treated"), "control": ("*", "control"), "controls": ("*", "control"),
    "untreated": ("treatment_role", "control"), "vehicle": ("treatment_role", "control"),
    "wild-type": ("genotype_role", "wild type"), "wildtype": ("genotype_role", "wild type"), "wt": ("genotype_role", "wild type"),
    "wild type": ("genotype_role", "wild type"), "knockout": ("genotype_role", "mutant"), "ko": ("genotype_role", "mutant"),
    "mutant": ("genotype_role", "mutant"), "mutants": ("genotype_role", "mutant"), "knockdown": ("genotype_role", "mutant"),
}
FILTER_STOP = {"the", "a", "an", "of", "which", "what", "are", "is", "samples", "sample", "files", "file", "in", "for",
               "specific", "and", "with", "from", "show", "me", "only", "all", "runs", "raw", "group"}


def project_summary(store: Store, accession: str) -> dict:
    data = store.load_projects([accession])
    if accession not in data:
        return {"error": f"{accession} is not in the index (only studies with curated SDRF are indexed)."}
    pd = data[accession]
    terms = store.sql_dicts("SELECT facet, value, n_samples FROM project_terms WHERE accession = ? ORDER BY facet, n_samples DESC", [accession])
    facets: dict[str, list] = defaultdict(list)
    for t in terms:
        facets[t["facet"]].append({"value": t["value"], "samples": t["n_samples"]})
    # experimental design: every factor / varying characteristic with its levels
    design = {}
    by_attr: dict[str, Counter] = defaultdict(Counter)
    for s, attrs in pd.samples.items():
        for a, vals in attrs.items():
            for v in vals:
                by_attr[a][v] += 1
    for a, cnt in sorted(by_attr.items()):
        if a in pd.factor_attrs or (1 < len(cnt) <= 30 and a not in vocab.INDIVIDUAL_ATTRS and a not in ("biological replicate",)):
            design[a] = {"factor": a in pd.factor_attrs, "levels": [{"value": v, "samples": n} for v, n in cnt.most_common(30)]}
    info = dict(pd.info)
    info.pop("doc_text", None)
    return {
        "accession": accession,
        "repository": info["repository"],
        "pride": _jsonable(pd.pride) if pd.pride else None,
        "counts": {
            "samples": info["n_samples"],
            "runs": info["n_runs"],
            "assay_rows": info["n_assays"],
            "individuals": info["n_individuals"],
            "individual_column": info["individual_column"],
            "max_fraction": info["max_fraction"],
            "label_channels": info["n_label_channels"],
        },
        "facets": dict(facets),
        "design": design,
        "modifications": [{k: m[k] for k in ("name", "unimod", "targets", "mod_type", "ptm_class")} for m in pd.modifications],
        "known_issues": pd.known_issues,
        "links": store.links(accession, info["sdrf_files"]),
        "provenance": {"sdrf_files": info["sdrf_files"], "sdrf_commit": store.sdrf_commit},
    }


def _jsonable(d):
    return {k: (str(v) if hasattr(v, "isoformat") else v) for k, v in d.items()}


def _has_word(text: str, phrase: str) -> bool:
    return re.search(r"(?<![a-z0-9])" + re.escape(phrase) + r"(?![a-z0-9])", text) is not None


def _sample_text(s) -> str:
    parts = [s.name] + [f"{a} {v}" for a, vals in s.attrs.items() for v in vals]
    return vocab.norm(" ".join(parts))


def sample_file_map(store: Store, accession: str, where: str | None = None, limit: int = 500) -> dict:
    """Samples (with attributes and roles) and their files, optionally filtered by a plain-language phrase."""
    data = store.load_projects([accession])
    if accession not in data:
        return {"error": f"{accession} is not in the index."}
    pd = data[accession]
    views = classify_samples(pd)
    notes = []
    selected = views
    if where and where.strip():
        w = " " + vocab.norm(where) + " "
        terms: list[tuple[str, object]] = []
        for phrase in sorted(ROLE_WORDS, key=len, reverse=True):
            if re.search(r"(?<![a-z0-9])" + re.escape(phrase) + r"(?![a-z0-9])", w):
                terms.append((phrase, ROLE_WORDS[phrase]))
                w = re.sub(r"(?<![a-z0-9])" + re.escape(phrase) + r"(?![a-z0-9])", " ", w)
        toks = [tok for tok in re.findall(r"[a-z0-9][a-z0-9\-\.+/]*", w) if tok not in FILTER_STOP]
        phrase = " ".join(toks)
        # prefer the whole remaining phrase ("patient 3"); fall back to individual words
        if len(toks) > 1 and any(_has_word(_sample_text(s), phrase) for s in selected):
            terms.append((phrase, None))
        else:
            terms.extend((tok, None) for tok in toks)
        for phrase, role in terms:
            if role is None:
                hit = [s for s in selected if _has_word(_sample_text(s), phrase)]
            else:
                kind, val = role
                hit = [s for s in selected if any(
                    r == val for k, (r, _) in getattr(s, "roles", {}).items() if kind == "*" or k == kind
                )]
            if hit:
                selected = hit
            else:
                notes.append(f"'{phrase}' matched no remaining samples, so it was ignored")
    runs = pd.runs_by_sample()
    rows = []
    for s in selected[:limit]:
        files = []
        for a in runs.get(s.name, []):
            url, src = _url_for(pd, a)
            files.append({
                "data_file": a["data_file"], "fraction": a["fraction"], "label": a["label"],
                "technical_replicate": a["technical_replicate"], "format": a["file_format"], "url": url, "url_source": src,
            })
        rows.append({
            "sample": s.name,
            "individual": s.individual,
            "roles": {k: {"role": r, "evidence": ev} for k, (r, ev) in getattr(s, "roles", {}).items()},
            "technical": s.technical,
            "cell_line": s.cell_line,
            "attributes": s.attrs,
            "files": files,
        })
    n_files = len({f["data_file"] for r in rows for f in r["files"]})
    return {
        "accession": accession,
        "filter": where,
        "notes": notes,
        "matched_samples": len(selected),
        "total_samples": len(views),
        "matched_files": n_files,
        "answer": f"{len(selected)} of {len(views)} samples match" + (f" '{where}'" if where else "") + f", mapped to {n_files} data files.",
        "samples": rows,
    }


def _url_for(pd, assay) -> tuple[str | None, str]:
    if assay.get("file_uri"):
        return assay["file_uri"], "SDRF comment[file uri]"
    pf = pd.pride_files.get(assay.get("data_file") or "")
    if pf and pf.get("url"):
        return pf["url"], "PRIDE file API"
    return None, "unresolved (PRIDE file API: /projects/{acc}/files/all)".replace("{acc}", pd.accession)

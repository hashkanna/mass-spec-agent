"""MCP server: exposes the cohort finder to Claude (or any MCP client) over stdio.

Add to Claude Code:     claude mcp add msagent -- msagent-mcp
Claude Desktop config:  {"mcpServers": {"msagent": {"command": "msagent-mcp"}}}
"""

from __future__ import annotations

import os
from functools import lru_cache

try:  # mcp >= 2
    from mcp.server.mcpserver import MCPServer as _Server
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP as _Server

from .cohort import find_cohorts as _find_cohorts
from .cohort import graph_query as _graph_query
from .criteria import Criteria, parse_query
from .manifest import build_manifest, to_tsv
from .project import project_summary, sample_file_map
from .store import Store

INSTRUCTIONS = """\
Mass-Spec Agent finds public proteomics datasets (PRIDE / ProteomeXchange) that can supply a cohort for
an experiment, using curated SDRF sample metadata (sample -> run -> file). Typical flow:
1. interpret_question to see the explicit criteria; edit them if they misread the question.
2. find_cohorts (or find_cohorts_structured with edited criteria). Every study comes back with each
   criterion marked confirmed / excluded / unknown plus the evidence. Report unknowns honestly;
   never present an unknown as confirmed.
3. list_files to answer "which files are the controls / tumours / WT samples" for one study.
4. export_manifest to hand the user a sample-to-file table with download URLs and provenance.
Counts distinguish individuals, biological samples, runs (data files) and fractions/label channels.
"""

mcp = _Server(name="msagent", instructions=INSTRUCTIONS)


@lru_cache(maxsize=1)
def store() -> Store:
    return Store()


def _compact(result: dict, max_checks_detail: int = 300) -> dict:
    for r in result.get("results", []):
        r.pop("cohort_preview", None)
        r.pop("facet_score", None)
        for ck in r.get("checks", []):
            ck["detail"] = ck["detail"][:max_checks_detail]
            for k in ("kind", "source"):
                ck.pop(k, None)
            if not ck.get("hint"):
                ck.pop("hint", None)
    return result


@mcp.tool()
def interpret_question(question: str) -> dict:
    """Turn a natural-language experimental question into explicit, editable search criteria.

    Returns the criteria object (pass it, edited if needed, to find_cohorts_structured) and a readable list.
    """
    c = parse_query(question, [v for v, _ in store().facet_values("organism")])
    return {"criteria": c.to_dict(), "criteria_list": c.describe()}


@mcp.tool()
def find_cohorts(question: str, limit: int = 10) -> dict:
    """Find public proteomics studies that satisfy an experimental question, verified criterion by criterion.

    Example question: "human liver tumour with paired non-tumour tissue, at least 10 patients, raw files,
    exclude cell lines". Results are ordered: qualifies, needs review (some criteria unknown), then near
    misses (excluded, with the reason).
    """
    return _compact(_find_cohorts(store(), question, limit=limit))


@mcp.tool()
def find_cohorts_structured(criteria: dict, limit: int = 10) -> dict:
    """Same as find_cohorts but with an explicit criteria object (as returned by interpret_question).

    Keys: organisms, tissues, diseases, cell_types, cell_lines, instruments, labels, tmt_plex, acquisition,
    ptms, comparison (tumour_vs_normal | disease_vs_control | treated_vs_control | mutant_vs_wildtype),
    paired, min_individuals, min_samples, exclude_cell_lines, only_cell_lines, raw_files, peak_lists,
    published_after, keywords, text.
    """
    return _compact(_find_cohorts(store(), Criteria.from_dict(criteria), limit=limit))


@mcp.tool()
def semantic_search(query: str, limit: int = 20) -> dict:
    """Hybrid keyword + embedding search over study metadata (no verification). Good for topics that are
    not structured fields, e.g. 'spliceosome', 'circadian', 'extracellular vesicles'."""
    st = store()
    if st.text is None:
        return {"error": "text index missing; run `msagent index`"}
    hits = st.text.hybrid(query, k=limit)
    out = []
    for acc, score in hits:
        _, n_samples, n_runs, year = st.projects[acc]
        out.append({"accession": acc, "score": round(score, 4), "samples": n_samples, "runs": n_runs})
    return {"query": query, "results": out}


@mcp.tool()
def graph_query(
    organisms: list[str] | None = None,
    tissues: list[str] | None = None,
    diseases: list[str] | None = None,
    cell_types: list[str] | None = None,
    instruments: list[str] | None = None,
    labels: list[str] | None = None,
    acquisition: list[str] | None = None,
    ptms: list[str] | None = None,
    min_samples: int | None = None,
    limit: int = 50,
) -> dict:
    """Exact facet filter over the project graph (no ranking by relevance). Values are normalised:
    organisms like 'homo sapiens'; tissues like 'plasma', 'liver', 'brain'; labels 'TMT'|'SILAC'|'iTRAQ'|'label-free';
    acquisition 'DIA'|'DDA'|'PRM'; ptms like 'phosphorylation', 'ubiquitination'; instruments like 'astral', 'timstof'."""
    return _graph_query(store(), organisms=organisms, tissues=tissues, diseases=diseases, cell_types=cell_types,
                        instruments=instruments, labels=labels, acquisition=acquisition, ptms=ptms,
                        min_samples=min_samples, limit=limit)


@mcp.tool()
def list_facet_values(facet: str, contains: str = "", limit: int = 50) -> dict:
    """Browse the controlled vocabulary actually used in the index, with project counts.
    Facets: organism, organism_part, disease, cell_type, cell_line, instrument, instrument_family, label_type,
    tmt_plex, acquisition, ptm, ptm_enrichment, enrichment, treatment, compound, genotype, file_format, factor."""
    vals = store().facet_values(facet)
    if contains:
        vals = [(v, n) for v, n in vals if contains.lower() in v.lower()]
    return {"facet": facet, "values": [{"value": v, "projects": n} for v, n in vals[:limit]]}


@mcp.tool()
def get_project(accession: str) -> dict:
    """Experimental design of one study: facets, factor levels, counts (individuals/samples/runs), modifications,
    known deposit issues, links to PRIDE and the SDRF."""
    return project_summary(store(), accession.strip().upper())


@mcp.tool()
def list_files(accession: str, where: str = "", limit: int = 200) -> dict:
    """Sample -> run -> file map for one study, optionally filtered by a phrase such as 'wild-type controls',
    'tumour', 'treated', 'patient 3' or 'fraction'. Each sample carries its inferred roles with evidence."""
    return sample_file_map(store(), accession.strip().upper(), where or None, limit=limit)


@mcp.tool()
def export_manifest(accessions: list[str], question: str = "", path: str = "") -> dict:
    """Build a cohort manifest (one row per sample x data file, with group, individual, URL and provenance)
    for the given studies, selecting samples according to the question. Writes a TSV if `path` is given."""
    man = build_manifest(store(), [a.strip().upper() for a in accessions], question)
    tsv = to_tsv(man)
    out = {"n_rows": man["n_rows"], "n_files": man["n_files"], "studies": [
        {k: s.get(k) for k in ("accession", "status", "design", "error")} for s in man["studies"]]}
    if path:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(tsv)
        out["path"] = os.path.abspath(path)
    else:
        out["tsv_preview"] = "\n".join(tsv.splitlines()[:40])
    return out


def main():
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()

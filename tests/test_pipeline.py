import pytest

from msagent.cohort import find_cohorts, graph_query
from msagent.manifest import build_manifest, to_tsv
from msagent.project import project_summary, sample_file_map


def by_acc(res):
    return {r["accession"]: r for r in res["results"]}


def check(r, cid):
    return next(c for c in r["checks"] if c["id"] == cid)


def test_ingest_counts(store):
    assert len(store.projects) == 4
    info = store.sql_dicts("SELECT * FROM projects WHERE accession = 'PXD990003'")[0]
    assert info["n_samples"] == 10 and info["n_runs"] == 2 and info["n_assays"] == 20
    assert info["n_label_channels"] == 10 and info["max_fraction"] == 2
    p1 = store.sql_dicts("SELECT * FROM projects WHERE accession = 'PXD990001'")[0]
    assert p1["n_individuals"] == 3 and p1["pub_year"] == 2021
    assert store.facets["tmt_plex"]["TMT 10/11-plex"].keys() == {"PXD990003"}


def test_paired_liver_cohort(store):
    res = find_cohorts(store, "human liver tumour with paired non-tumour tissue, at least 3 patients, raw files, exclude cell lines")
    r = by_acc(res)
    assert r["PXD990001"]["status"] == "qualifies"
    assert check(r["PXD990001"], "paired")["status"] == "confirmed"
    assert r["PXD990001"]["counts"]["individuals"] == 3
    assert r["PXD990001"]["counts"]["groups"] == {"tumour": 3, "normal": 3}
    # tumour-only study: the comparison is excluded with a reason, and it is flagged by curators
    other = r["PXD990002"]
    assert other["status"] == "excluded"
    assert check(other, "comparison")["status"] == "excluded"
    assert check(other, "known_issues")["status"] == "excluded"
    assert res["results"][0]["accession"] == "PXD990001"


def test_individual_threshold_is_enforced(store):
    res = find_cohorts(store, "human liver tumour with paired normal tissue, at least 5 patients")
    r = by_acc(res)["PXD990001"]
    assert r["status"] == "excluded"
    assert check(r, "min_individuals")["status"] == "excluded"


def test_cell_line_ptm_treatment(store):
    res = find_cohorts(store, "TMT phosphoproteomics of HeLa cells treated with trametinib versus DMSO")
    r = by_acc(res)["PXD990003"]
    assert r["status"] == "qualifies", r["checks"]
    assert check(r, "ptm")["status"] == "confirmed"
    assert check(r, "comparison")["detail"].startswith("5 treated")
    # the same study is excluded when cell lines are not wanted
    res = find_cohorts(store, "human phosphoproteomics, exclude cell lines")
    assert check(by_acc(res)["PXD990003"], "no_cell_lines")["status"] == "excluded"


def test_unknowns_are_not_confirmed(store):
    res = find_cohorts(store, "mouse brain DIA timsTOF knockout vs wild type, at least 3 individuals, published since 2024")
    r = by_acc(res)["PXD990004"]
    assert check(r, "comparison")["status"] == "confirmed"
    assert check(r, "min_individuals")["status"] == "unknown"  # no individual column, 6 samples
    assert check(r, "published_after")["status"] == "unknown"  # no PRIDE date
    assert r["status"] == "needs review"


def test_sample_file_map_filters(store):
    m = sample_file_map(store, "PXD990004", "wild-type")
    assert m["matched_samples"] == 3
    assert all(s["roles"]["genotype_role"]["role"] == "wild type" for s in m["samples"])
    m = sample_file_map(store, "PXD990001", "patient 2 tumour")
    assert [s["sample"] for s in m["samples"]] == ["patient2_T"]
    assert m["samples"][0]["files"][0]["url"].endswith("/PXD990001/P2_T.raw")


def test_manifest_and_summary(store):
    man = build_manifest(store, ["PXD990001"], "liver tumour with paired normal tissue")
    assert man["n_rows"] == 6 and man["n_files"] == 6
    assert {r["group"] for r in man["rows"]} == {"tumour", "normal"}
    tsv = to_tsv(man)
    assert "PXD990001\ttumour" in tsv and tsv.startswith("# msagent cohort manifest")
    s = project_summary(store, "PXD990003")
    assert "treatment" in s["design"] and s["counts"]["label_channels"] == 10


def test_graph_query(store):
    out = graph_query(store, organisms=["homo sapiens"], tissues=["liver"])
    assert {r["accession"] for r in out["results"]} == {"PXD990001", "PXD990002"}


def test_mcp_tools_registered():
    pytest.importorskip("mcp")
    import asyncio

    from msagent.mcp_server import mcp

    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert {"find_cohorts", "semantic_search", "graph_query", "list_files", "export_manifest"} <= names

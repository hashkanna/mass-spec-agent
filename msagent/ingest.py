"""Build the DuckDB metadata graph from curated SDRF files.

Grain of the tables:

* ``projects``           one row per accession (PXD/MSV/...), with counts and provenance
* ``assays``             one row per SDRF row (sample x data file x label channel x fraction)
* ``sample_attributes``  long table of sample characteristics and factor values
* ``project_terms``      normalised facet values per project with sample counts (the filter index)
* ``modifications``      searched modifications per project
* ``known_issues``       deposits flagged by curators (docs/known-issues.md)
* ``pride_projects`` / ``pride_files``  optional enrichment from the PRIDE Archive API

Spectra are never touched; only metadata and file locations are stored.
"""

from __future__ import annotations

import re
import subprocess
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pyarrow as pa

from . import vocab
from .config import DB_PATH, SDRF_DIR
from .sdrf import is_missing, iter_sdrf_paths, parse_cell, read_sdrf

SCHEMA = """
CREATE TABLE projects (
    accession VARCHAR PRIMARY KEY,
    repository VARCHAR,
    sdrf_files VARCHAR[],
    n_assays INTEGER,
    n_samples INTEGER,
    n_runs INTEGER,
    n_individuals INTEGER,
    individual_column VARCHAR,
    max_fraction INTEGER,
    n_label_channels INTEGER,
    templates VARCHAR[],
    file_uri_rows INTEGER,
    pub_year INTEGER,
    pub_month INTEGER,
    doc_text VARCHAR
);
CREATE TABLE assays (
    accession VARCHAR,
    sdrf_file VARCHAR,
    row_id INTEGER,
    source_name VARCHAR,
    assay_name VARCHAR,
    data_file VARCHAR,
    file_uri VARCHAR,
    file_format VARCHAR,
    fraction VARCHAR,
    technical_replicate VARCHAR,
    biological_replicate VARCHAR,
    label VARCHAR,
    instrument VARCHAR,
    acquisition VARCHAR
);
CREATE TABLE sample_attributes (
    accession VARCHAR,
    source_name VARCHAR,
    kind VARCHAR,
    attribute VARCHAR,
    value VARCHAR,
    term_accession VARCHAR
);
CREATE TABLE project_terms (
    accession VARCHAR,
    facet VARCHAR,
    value VARCHAR,
    n_samples INTEGER,
    n_assays INTEGER
);
CREATE TABLE modifications (
    accession VARCHAR,
    name VARCHAR,
    unimod VARCHAR,
    targets VARCHAR,
    mod_type VARCHAR,
    position VARCHAR,
    ptm_class VARCHAR
);
CREATE TABLE known_issues (
    accession VARCHAR,
    severity VARCHAR,
    category VARCHAR,
    section VARCHAR,
    detail VARCHAR
);
CREATE TABLE meta (key VARCHAR PRIMARY KEY, value VARCHAR);
CREATE TABLE pride_projects (
    accession VARCHAR PRIMARY KEY,
    title VARCHAR,
    description VARCHAR,
    sample_protocol VARCHAR,
    data_protocol VARCHAR,
    keywords VARCHAR[],
    submission_date DATE,
    publication_date DATE,
    pubmed_ids VARCHAR[],
    dois VARCHAR[],
    raw_json VARCHAR
);
CREATE TABLE pride_files (
    accession VARCHAR,
    file_name VARCHAR,
    category VARCHAR,
    size_bytes BIGINT,
    url VARCHAR,
    checksum VARCHAR
);
"""

# sample characteristic -> facet name
SAMPLE_FACETS = {
    "organism": "organism",
    "disease": "disease",
    "cell type": "cell_type",
    "cell line": "cell_line",
    "cellosaurus name": "cell_line",
    "sex": "sex",
    "developmental stage": "developmental_stage",
    "material type": "material_type",
    "sample type": "sample_type",
    "enrichment process": "enrichment",
    "treatment": "treatment",
    "compound": "compound",
    "genotype": "genotype",
    "phenotype": "phenotype",
    "strain or breed": "strain",
}
for _a in vocab.TISSUE_ATTRS:
    SAMPLE_FACETS[_a] = "organism_part"

_FTP_DATE = re.compile(r"/archive/(\d{4})/(\d{2})/")


def repository_of(accession: str) -> str:
    a = accession.upper()
    for prefix, repo in [
        ("PXD", "ProteomeXchange/PRIDE"),
        ("PAD", "PRIDE (affinity)"),
        ("MSV", "MassIVE"),
        ("IPX", "iProX"),
        ("JPST", "jPOST"),
        ("PDC", "Proteomic Data Commons"),
        ("PRD", "PRIDE"),
    ]:
        if a.startswith(prefix):
            return repo
    return "other"


def _ival(s: str | None) -> int | None:
    try:
        return int(float(s)) if s is not None else None
    except ValueError:
        return None


class _ProjectAcc:
    """Accumulates everything for one accession across its SDRF files."""

    def __init__(self, accession: str):
        self.accession = accession
        self.sdrf_files: list[str] = []
        self.templates: set[str] = set()
        self.assays: list[tuple] = []
        self.attrs: set[tuple] = set()
        self.term_samples: dict[tuple[str, str], set[str]] = defaultdict(set)
        self.term_assays: Counter = Counter()
        self.mods: dict[tuple, vocab.Modification] = {}
        self.samples: set[str] = set()
        self.runs: set[str] = set()
        self.individual_cols: Counter = Counter()
        self.individuals: set[str] = set()
        self.max_fraction = 0
        self.labels: set[str] = set()
        self.file_uri_rows = 0
        self.pub: Counter = Counter()


def _add_file(acc: _ProjectAcc, rel_path: str, sdrf) -> None:
    acc.sdrf_files.append(rel_path)
    row_offset = len(acc.assays)
    for i, r in enumerate(sdrf.rows):
        src = r.source_name or r.assay_name or f"row{i + 1}"
        fmt = vocab.ext_format(r.data_file) if r.data_file else None
        acc.assays.append(
            (
                acc.accession, rel_path, row_offset + i + 1, src, r.assay_name, r.data_file, r.file_uri, fmt,
                r.fraction, r.technical_replicate, r.biological_replicate, r.label, r.instrument, r.acquisition,
            )
        )
        acc.samples.add(src)
        if r.data_file:
            acc.runs.add(r.data_file)
        if r.file_uri:
            acc.file_uri_rows += 1
            m = _FTP_DATE.search(r.file_uri)
            if m:
                acc.pub[(int(m.group(1)), int(m.group(2)))] += 1
        f = _ival(r.fraction)
        # some SDRFs put identifiers rather than indices here; only small integers are fraction numbers
        if f is not None and acc.max_fraction < f <= 1000:
            acc.max_fraction = f

        # sample-level attributes
        for kind, block in (("characteristics", r.characteristics), ("factor", r.factors)):
            for attr, values in block.items():
                for value, term_acc in values:
                    if is_missing(value):
                        continue
                    acc.attrs.add((acc.accession, src, kind, attr, value, term_acc))
                    if kind == "characteristics" and attr in vocab.INDIVIDUAL_ATTRS:
                        acc.individual_cols[attr] += 1
                        if vocab.norm(value) not in {"pooled", "pool", "mixed"}:
                            acc.individuals.add(f"{attr}={value}")
                    facet = SAMPLE_FACETS.get(attr)
                    if facet is None:
                        continue
                    if facet == "organism":
                        val = vocab.normalize_organism(value)
                    else:
                        val = vocab.norm(value)
                    acc.term_samples[(facet, val)].add(src)
                    if facet == "enrichment":
                        ptm = vocab.enrichment_ptm(value)
                        if ptm:
                            acc.term_samples[("ptm_enrichment", ptm)].add(src)
            if kind == "factor":
                for attr in block:
                    acc.term_samples[("factor", attr)].add(src)

        # run-level technical facets
        tech: list[tuple[str, str]] = []
        inst = vocab.instrument_info(r.instrument)
        if inst:
            tech += [("instrument", inst[0]), ("instrument_family", inst[1]), ("vendor", inst[2])]
        lt = vocab.label_type(r.label)
        if lt:
            tech.append(("label_type", lt))
            acc.labels.add(r.label)
        aq = vocab.acquisition_type(r.acquisition)
        if aq:
            tech.append(("acquisition", aq))
        if fmt:
            tech.append(("file_format", fmt))
        for cleave in r.comments.get("cleavage agent details", []):
            if not is_missing(cleave):
                tech.append(("cleavage_agent", parse_cell(cleave)["NT"]))
        for frac in r.comments.get("fractionation method", []):
            if not is_missing(frac):
                tech.append(("fractionation", vocab.norm(parse_cell(frac)["NT"])))
        for tmpl in r.comments.get("sdrf template", []):
            if not is_missing(tmpl):
                acc.templates.add(parse_cell(tmpl)["NT"])
        for mp in r.comments.get("modification parameters", []):
            if is_missing(mp):
                continue
            mod = vocab.parse_modification(parse_cell(mp))
            if mod:
                key = (vocab.norm(mod.name), mod.targets, mod.mod_type, mod.position)
                acc.mods.setdefault(key, mod)
        for facet, val in tech:
            acc.term_samples[(facet, val)].add(src)
            acc.term_assays[(facet, val)] += 1


def _finalise(acc: _ProjectAcc):
    terms = {k: v for k, v in acc.term_samples.items()}
    # project-level derived facets
    plex = vocab.tmt_plex(acc.labels)
    mod_names = {vocab.norm(m.name) for m in acc.mods.values()}
    if "tmtpro" in mod_names or "tmtpro16plex" in mod_names:
        plex = "TMTpro (16/18-plex)"
    if plex:
        terms[("tmt_plex", plex)] = set(terms.get(("label_type", "TMT"), set()))
    for m in acc.mods.values():
        cls = m.ptm_class
        terms.setdefault(("modification", m.name), set())
        if cls:
            terms.setdefault(("ptm", cls), set())
    if acc.labels and not any(k[0] == "label_type" and k[1] == "TMT" for k in terms):
        if any(n.startswith("tmt") for n in mod_names):
            terms[("label_type", "TMT")] = set()
    for tmpl in acc.templates:
        terms[("template", tmpl)] = set()

    term_rows = [
        (acc.accession, facet, val, len(samples), acc.term_assays.get((facet, val), 0))
        for (facet, val), samples in terms.items()
    ]
    ind_col = acc.individual_cols.most_common(1)[0][0] if acc.individual_cols else None
    pub = acc.pub.most_common(1)[0][0] if acc.pub else (None, None)
    project_row = dict(
        accession=acc.accession,
        repository=repository_of(acc.accession),
        sdrf_files=acc.sdrf_files,
        n_assays=len(acc.assays),
        n_samples=len(acc.samples),
        n_runs=len(acc.runs),
        n_individuals=len(acc.individuals) if ind_col else None,
        individual_column=ind_col,
        max_fraction=acc.max_fraction or None,
        n_label_channels=len(acc.labels) or None,
        templates=sorted(acc.templates),
        file_uri_rows=acc.file_uri_rows,
        pub_year=pub[0],
        pub_month=pub[1],
        doc_text=build_doc_text(acc.accession, terms, acc),
    )
    mod_rows = [
        (acc.accession, m.name, m.accession, m.targets, m.mod_type, m.position, m.ptm_class) for m in acc.mods.values()
    ]
    return project_row, term_rows, mod_rows


_DOC_FACETS = [
    ("organism", "Organism"),
    ("organism_part", "Organism part / tissue"),
    ("disease", "Disease"),
    ("cell_type", "Cell type"),
    ("cell_line", "Cell line"),
    ("material_type", "Material"),
    ("developmental_stage", "Developmental stage"),
    ("treatment", "Treatment"),
    ("compound", "Compound"),
    ("genotype", "Genotype"),
    ("phenotype", "Phenotype"),
    ("enrichment", "Enrichment"),
    ("ptm", "PTMs searched"),
    ("instrument", "Instrument"),
    ("label_type", "Quantification"),
    ("tmt_plex", "TMT plex"),
    ("acquisition", "Acquisition"),
    ("fractionation", "Fractionation"),
    ("sample_type", "Sample type"),
    ("template", "SDRF template"),
    ("factor", "Experimental factors"),
]


def build_doc_text(accession: str, terms: dict, acc: _ProjectAcc) -> str:
    by_facet: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for (facet, val), samples in terms.items():
        by_facet[facet].append((val, len(samples)))
    parts = [f"{accession}."]
    for facet, label in _DOC_FACETS:
        vals = sorted(by_facet.get(facet, []), key=lambda x: -x[1])[:12]
        if vals:
            parts.append(f"{label}: " + "; ".join(v for v, _ in vals) + ".")
    mods = sorted({m.name for m in acc.mods.values()})
    if mods:
        parts.append("Modifications: " + ", ".join(mods[:15]) + ".")
    # sample names often encode the design (e.g. "tumor_patient3", "WT_ctrl_rep1")
    names = sorted(acc.samples)[:25]
    parts.append("Samples: " + ", ".join(names) + ".")
    return " ".join(parts)


def parse_known_issues(path: Path) -> list[tuple]:
    if not path.exists():
        return []
    rows = []
    severity = section = category = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"^##\s+(Critical|Major|Moderate|Minor)\s*:\s*(.*)$", line.strip(), re.I)
        if m:
            severity, section, category = m.group(1).lower(), m.group(2).strip(), None
            continue
        m = re.match(r"^Category:\s*`([^`]+)`", line.strip())
        if m:
            category = m.group(1)
            continue
        m = re.match(r"^\|\s*((?:PXD|MSV|IPX|PAD|JPST)\d+)\s*\|(.*)$", line.strip())
        if m and severity:
            detail = " | ".join(c.strip() for c in m.group(2).split("|") if c.strip())
            rows.append((m.group(1), severity, category, section, detail[:600]))
    return rows


def _git_commit(path: Path) -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(path), "log", "-1", "--format=%H %cI"], capture_output=True, text=True, timeout=10)
        return out.stdout.strip() or None
    except Exception:
        return None


def build(sdrf_dir: Path = SDRF_DIR, db_path: Path = DB_PATH, limit: int | None = None, verbose: bool = True) -> dict:
    sdrf_dir = Path(sdrf_dir)
    if not (sdrf_dir / "datasets").exists():
        raise SystemExit(
            f"SDRF corpus not found at {sdrf_dir}. Run scripts/fetch_data.sh (or set MSAGENT_SDRF_DIR)."
        )
    t0 = time.time()
    paths = list(iter_sdrf_paths(sdrf_dir))
    by_acc: dict[str, list[Path]] = defaultdict(list)
    for acc, p in paths:
        by_acc[acc].append(p)
    accessions = sorted(by_acc)
    if limit:
        accessions = accessions[:limit]

    projects, terms, mods = [], [], []
    assay_batches: list[list[tuple]] = []
    attr_batches: list[list[tuple]] = []
    failures = []
    for i, accession in enumerate(accessions):
        acc = _ProjectAcc(accession)
        for p in by_acc[accession]:
            try:
                _add_file(acc, str(p.relative_to(sdrf_dir)), read_sdrf(p, accession))
            except Exception as exc:  # keep going; report at the end
                failures.append((str(p), repr(exc)))
        if not acc.assays:
            continue
        prow, trows, mrows = _finalise(acc)
        projects.append(prow)
        terms.extend(trows)
        mods.extend(mrows)
        assay_batches.append(acc.assays)
        attr_batches.append(sorted(acc.attrs, key=lambda x: (x[1], x[2], x[3], x[4])))
        if verbose and (i + 1) % 1000 == 0:
            print(f"  parsed {i + 1}/{len(accessions)} accessions ({time.time() - t0:.0f}s)", flush=True)

    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = db_path.with_suffix(".building")
    if tmp.exists():
        tmp.unlink()
    con = duckdb.connect(str(tmp))
    con.execute(SCHEMA)

    def _insert(table: str, columns: list[str], rows_iter):
        cols = {c: [] for c in columns}
        for row in rows_iter:
            for c, v in zip(columns, row):
                cols[c].append(v)
        arrow = pa.table(cols)  # noqa: F841  (referenced by name inside SQL)
        con.execute(f"INSERT INTO {table} SELECT * FROM arrow")

    _insert(
        "projects",
        list(projects[0].keys()),
        (tuple(p.values()) for p in projects),
    )
    _insert(
        "assays",
        ["accession", "sdrf_file", "row_id", "source_name", "assay_name", "data_file", "file_uri", "file_format",
         "fraction", "technical_replicate", "biological_replicate", "label", "instrument", "acquisition"],
        (r for batch in assay_batches for r in batch),
    )
    _insert(
        "sample_attributes",
        ["accession", "source_name", "kind", "attribute", "value", "term_accession"],
        (r for batch in attr_batches for r in batch),
    )
    _insert("project_terms", ["accession", "facet", "value", "n_samples", "n_assays"], terms)
    _insert("modifications", ["accession", "name", "unimod", "targets", "mod_type", "position", "ptm_class"], mods)
    ki = parse_known_issues(sdrf_dir / "docs" / "known-issues.md")
    if ki:
        _insert("known_issues", ["accession", "severity", "category", "section", "detail"], ki)

    meta = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sdrf_source": "https://github.com/bigbio/sdrf-annotated-datasets",
        "sdrf_commit": _git_commit(sdrf_dir) or "unknown",
        "n_projects": str(len(projects)),
        "parse_failures": str(len(failures)),
    }
    con.executemany("INSERT INTO meta VALUES (?, ?)", list(meta.items()))
    con.execute("CREATE INDEX idx_terms ON project_terms(facet, value)")
    con.execute("CREATE INDEX idx_attr_acc ON sample_attributes(accession)")
    con.execute("CREATE INDEX idx_assay_acc ON assays(accession)")
    con.close()
    if db_path.exists():
        db_path.unlink()
    tmp.rename(db_path)
    stats = {
        "projects": len(projects),
        "assays": sum(len(b) for b in assay_batches),
        "sample_attributes": sum(len(b) for b in attr_batches),
        "terms": len(terms),
        "failures": failures[:20],
        "n_failures": len(failures),
        "seconds": round(time.time() - t0, 1),
        "db": str(db_path),
    }
    return stats

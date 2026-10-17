"""Read-side access to the metadata graph: facet index in memory, per-project loaders."""

from __future__ import annotations

import threading
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import duckdb

from .config import DB_PATH, INDEX_DIR, PRIDE_PROJECT_PAGE
from .index import TextIndex

SDRF_BLOB = "https://github.com/bigbio/sdrf-annotated-datasets/blob/{commit}/{path}"


@dataclass
class ProjectData:
    accession: str
    info: dict
    samples: dict[str, dict[str, list[str]]]  # source name -> attribute -> values (characteristics + factor values)
    factor_attrs: set[str]
    assays: list[dict]
    modifications: list[dict]
    known_issues: list[dict]
    pride: dict | None = None
    pride_files: dict[str, dict] = field(default_factory=dict)

    def runs_by_sample(self) -> dict[str, list[dict]]:
        out: dict[str, list[dict]] = defaultdict(list)
        for a in self.assays:
            out[a["source_name"]].append(a)
        return out


class Store:
    def __init__(self, db_path: Path = DB_PATH, index_dir: Path = INDEX_DIR, load_text_index: bool = True):
        if not Path(db_path).exists():
            raise FileNotFoundError(f"No database at {db_path}; run `msagent ingest` first.")
        self._con = duckdb.connect(str(db_path), read_only=True)
        self._lock = threading.Lock()
        self.meta = dict(self.sql("SELECT key, value FROM meta"))
        commit = (self.meta.get("sdrf_commit") or "main").split()[0]
        self.sdrf_commit = commit
        # facet -> value -> {accession: n_samples}
        self.facets: dict[str, dict[str, dict[str, int]]] = defaultdict(lambda: defaultdict(dict))
        for acc, facet, value, n in self.sql("SELECT accession, facet, value, n_samples FROM project_terms"):
            self.facets[facet][value][acc] = n
        self.annotated: dict[str, set[str]] = {
            f: {a for accs in vals.values() for a in accs} for f, vals in self.facets.items()
        }
        self.projects = {r[0]: r for r in self.sql("SELECT accession, n_samples, n_runs, pub_year FROM projects")}
        self.pride_years = dict(self.sql("SELECT accession, year(publication_date) FROM pride_projects WHERE publication_date IS NOT NULL"))
        self.text: TextIndex | None = None
        if load_text_index:
            try:
                self.text = TextIndex(index_dir)
            except FileNotFoundError:
                self.text = None

    # duckdb connections are not safe to share across threads; serialise access
    def sql(self, query: str, params: list | None = None) -> list[tuple]:
        with self._lock:
            return self._con.execute(query, params or []).fetchall()

    def sql_dicts(self, query: str, params: list | None = None) -> list[dict]:
        with self._lock:
            cur = self._con.execute(query, params or [])
            cols = [d[0] for d in cur.description]
            return [dict(zip(cols, r)) for r in cur.fetchall()]

    def facet_values(self, facet: str) -> list[tuple[str, int]]:
        """Distinct values of a facet with the number of projects using them."""
        vals = self.facets.get(facet, {})
        return sorted(((v, len(a)) for v, a in vals.items()), key=lambda x: -x[1])

    def sdrf_url(self, rel_path: str) -> str:
        return SDRF_BLOB.format(commit=self.sdrf_commit, path=rel_path)

    def links(self, accession: str, sdrf_files: list[str]) -> dict:
        out = {"sdrf": [self.sdrf_url(p) for p in sdrf_files]}
        if accession.startswith(("PXD", "PAD")):
            out["pride"] = PRIDE_PROJECT_PAGE.format(accession=accession)
        if accession.startswith("PXD"):
            out["proteomexchange"] = f"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={accession}"
        if accession.startswith("MSV"):
            out["massive"] = f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession={accession}"
        return out

    def load_projects(self, accessions: list[str]) -> dict[str, ProjectData]:
        if not accessions:
            return {}
        ph = ",".join("?" for _ in accessions)
        info = {r["accession"]: r for r in self.sql_dicts(f"SELECT * FROM projects WHERE accession IN ({ph})", accessions)}
        samples: dict[str, dict[str, dict[str, list[str]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        factor_attrs: dict[str, set[str]] = defaultdict(set)
        for acc, src, kind, attr, value in self.sql(
            f"SELECT accession, source_name, kind, attribute, value FROM sample_attributes WHERE accession IN ({ph})",
            accessions,
        ):
            vals = samples[acc][src][attr]
            if value not in vals:
                vals.append(value)
            if kind == "factor":
                factor_attrs[acc].add(attr)
        assays: dict[str, list[dict]] = defaultdict(list)
        for row in self.sql_dicts(f"SELECT * FROM assays WHERE accession IN ({ph}) ORDER BY accession, row_id", accessions):
            assays[row["accession"]].append(row)
        mods: dict[str, list[dict]] = defaultdict(list)
        for row in self.sql_dicts(f"SELECT * FROM modifications WHERE accession IN ({ph})", accessions):
            mods[row["accession"]].append(row)
        issues: dict[str, list[dict]] = defaultdict(list)
        for row in self.sql_dicts(f"SELECT * FROM known_issues WHERE accession IN ({ph})", accessions):
            issues[row["accession"]].append(row)
        pride = {r["accession"]: r for r in self.sql_dicts(
            f"SELECT accession, title, description, sample_protocol, data_protocol, keywords, submission_date, publication_date, pubmed_ids, dois FROM pride_projects WHERE accession IN ({ph})",
            accessions,
        )}
        pfiles: dict[str, dict[str, dict]] = defaultdict(dict)
        for row in self.sql_dicts(f"SELECT * FROM pride_files WHERE accession IN ({ph})", accessions):
            pfiles[row["accession"]][row["file_name"]] = row
        out = {}
        for acc in accessions:
            if acc not in info:
                continue
            # make sure samples that only appear in assays still exist
            for a in assays[acc]:
                samples[acc][a["source_name"]]
            out[acc] = ProjectData(
                accession=acc,
                info=info[acc],
                samples={s: dict(v) for s, v in samples[acc].items()},
                factor_attrs=factor_attrs[acc],
                assays=assays[acc],
                modifications=mods[acc],
                known_issues=issues[acc],
                pride=pride.get(acc),
                pride_files=pfiles.get(acc, {}),
            )
        return out

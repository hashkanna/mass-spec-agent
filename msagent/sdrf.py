"""Parsing of SDRF-Proteomics files (sample-to-data-file relationship tables).

An SDRF row is one *assay*: one sample (``source name``) measured in one data file
(``comment[data file]``), possibly one fraction / label channel / technical replicate of it.
We keep that grain, and split columns into sample-level characteristics, factor values and
run-level technical comments.
"""

from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

csv.field_size_limit(sys.maxsize)

_HEADER_RE = re.compile(r"^\s*(characteristics|comment|factor value|factorvalue)\s*\[\s*(.*?)\s*\]\s*$", re.I)
_KV_RE = re.compile(r"(?:^|;)\s*(NT|AC|TA|MT|PP|MM|CS|TS)\s*=", re.I)

MISSING = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "not available",
    "not applicable",
    "not apllicable",
    "not specified",
    "unknown",
    "not reported",
    "missing",
    "-",
    "--",
}


def is_missing(value: str | None) -> bool:
    return value is None or value.strip().lower() in MISSING


def parse_header(col: str) -> tuple[str, str]:
    """Return (kind, name) where kind is characteristics/comment/factor/core."""
    m = _HEADER_RE.match(col)
    if not m:
        return "core", " ".join(col.strip().lower().split())
    kind = m.group(1).lower().replace(" ", "")
    kind = {"characteristics": "characteristics", "comment": "comment"}.get(kind, "factor")
    return kind, " ".join(m.group(2).strip().lower().split())


def parse_cell(value: str) -> dict[str, str]:
    """Parse ``NT=Phospho;AC=UNIMOD:21;TA=S,T`` style cells.

    Plain cells come back as ``{"NT": value}``. Keys are upper-cased.
    """
    value = value.strip()
    if not _KV_RE.search(value):
        return {"NT": value}
    out: dict[str, str] = {}
    for part in value.split(";"):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        k = k.strip().upper()
        if k in out:  # e.g. TA listed twice
            out[k] = out[k] + "," + v.strip()
        else:
            out[k] = v.strip()
    if "NT" not in out:
        out["NT"] = out.get("AC", value)
    return out


def cell_name(value: str) -> str:
    return parse_cell(value)["NT"].strip()


@dataclass
class AssayRow:
    source_name: str
    assay_name: str
    data_file: str
    file_uri: str | None
    fraction: str | None
    technical_replicate: str | None
    biological_replicate: str | None
    label: str | None
    instrument: str | None
    acquisition: str | None
    characteristics: dict[str, list[tuple[str, str | None]]] = field(default_factory=dict)
    factors: dict[str, list[tuple[str, str | None]]] = field(default_factory=dict)
    comments: dict[str, list[str]] = field(default_factory=dict)


@dataclass
class SdrfFile:
    path: Path
    accession: str
    columns: list[tuple[str, str]]
    rows: list[AssayRow]


def _first(d: dict[str, list], key: str) -> str | None:
    vals = d.get(key) or []
    for v in vals:
        v = v[0] if isinstance(v, tuple) else v
        if not is_missing(v):
            return v
    return None


def read_sdrf(path: Path, accession: str | None = None) -> SdrfFile:
    path = Path(path)
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        try:
            header = next(reader)
        except StopIteration:
            return SdrfFile(path, accession or path.parent.name, [], [])
        columns = [parse_header(h) for h in header]
        rows: list[AssayRow] = []
        for raw in reader:
            if not raw or all(not c.strip() for c in raw):
                continue
            chars: dict[str, list[tuple[str, str | None]]] = {}
            factors: dict[str, list[tuple[str, str | None]]] = {}
            comments: dict[str, list[str]] = {}
            core: dict[str, str] = {}
            for (kind, name), val in zip(columns, raw):
                val = val.strip()
                if kind == "core":
                    core.setdefault(name, val)
                elif kind == "comment":
                    comments.setdefault(name, []).append(val)
                else:
                    parsed = parse_cell(val)
                    target = chars if kind == "characteristics" else factors
                    target.setdefault(name, []).append((parsed["NT"], parsed.get("AC")))
            label = _first(comments, "label")
            rows.append(
                AssayRow(
                    source_name=core.get("source name", ""),
                    assay_name=core.get("assay name", ""),
                    data_file=_first(comments, "data file") or "",
                    file_uri=_first(comments, "file uri"),
                    fraction=_first(comments, "fraction identifier"),
                    technical_replicate=_first(comments, "technical replicate"),
                    biological_replicate=_first(chars, "biological replicate"),
                    label=cell_name(label) if label else None,
                    instrument=cell_name(_first(comments, "instrument") or "") or None,
                    acquisition=cell_name(_first(comments, "proteomics data acquisition method") or "") or None,
                    characteristics=chars,
                    factors=factors,
                    comments=comments,
                )
            )
    return SdrfFile(path, accession or path.parent.name, columns, rows)


def iter_sdrf_paths(root: Path):
    """Yield (accession, path) for every curated SDRF under ``datasets/``."""
    ds = Path(root) / "datasets"
    for d in sorted(ds.iterdir()):
        if d.is_dir():
            for f in sorted(d.iterdir()):
                if f.is_file() and ".sdrf" in f.name:
                    yield d.name, f
        elif d.is_file() and ".sdrf" in d.name:
            yield d.name.split(".")[0].split("_")[0], d

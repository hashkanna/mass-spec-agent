"""2-D map of studies (UMAP over the project embeddings), coloured by metadata in the web UI.

The embeddings come from SDRF-derived metadata text (plus PRIDE text when enriched), so the map shows
which studies look alike in their annotated design, not in their spectra.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import duckdb
import numpy as np

from .config import DB_PATH, INDEX_DIR

COLOR_FACETS = ["organism", "instrument_family", "label_type", "acquisition", "organism_part", "disease", "ptm"]


def build(db_path: Path = DB_PATH, index_dir: Path = INDEX_DIR) -> dict:
    import umap

    index_dir = Path(index_dir)
    accs = json.loads((index_dir / "accessions.json").read_text())
    vecs = np.load(index_dir / "embeddings.npy")
    xy = umap.UMAP(n_neighbors=20, min_dist=0.15, metric="cosine", random_state=7).fit_transform(vecs)
    con = duckdb.connect(str(db_path), read_only=True)
    top: dict[str, dict[str, str]] = defaultdict(dict)
    for acc, facet, value in con.execute(
        f"""SELECT accession, facet, arg_max(value, n_samples) FROM project_terms
            WHERE facet IN ({','.join('?' for _ in COLOR_FACETS)}) GROUP BY accession, facet""",
        COLOR_FACETS,
    ).fetchall():
        top[acc][facet] = value
    sizes = dict(con.execute("SELECT accession, n_samples FROM projects").fetchall())
    con.close()
    xy = (xy - xy.min(0)) / np.maximum(xy.max(0) - xy.min(0), 1e-9)
    points = [
        {"a": a, "x": round(float(x), 4), "y": round(float(y), 4), "n": sizes.get(a, 0), **{f: top[a].get(f) for f in COLOR_FACETS}}
        for a, (x, y) in zip(accs, xy)
    ]
    out = {"facets": COLOR_FACETS, "points": points}
    (index_dir / "map.json").write_text(json.dumps(out, separators=(",", ":")))
    return {"points": len(points), "path": str(index_dir / "map.json")}

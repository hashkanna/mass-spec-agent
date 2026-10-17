"""Web UI + JSON API. Run with `msagent serve`."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from fastapi import Body, FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse

from ..cohort import find_cohorts
from ..config import INDEX_DIR, ROOT
from ..criteria import Criteria, parse_query
from ..manifest import build_manifest, to_json, to_tsv
from ..project import project_summary, sample_file_map
from ..store import Store

STATIC = Path(__file__).parent / "static"
app = FastAPI(title="Mass-Spec Agent", version="0.1.0")


@lru_cache(maxsize=1)
def store() -> Store:
    return Store()


def _json(obj) -> JSONResponse:
    return JSONResponse(json.loads(json.dumps(obj, default=str)))


@app.get("/")
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/api/status")
def status():
    st = store()
    return _json({"meta": st.meta, "projects": len(st.projects), "text_index": st.text.meta if st.text else None})


@app.get("/api/criteria")
def criteria(q: str):
    c = parse_query(q, [v for v, _ in store().facet_values("organism")])
    return _json({"criteria": c.to_dict(), "criteria_list": c.describe()})


@app.get("/api/search")
def search(q: str, limit: int = 15):
    return _json(find_cohorts(store(), q, limit=limit))


@app.post("/api/search")
def search_structured(criteria: dict = Body(...), limit: int = 15):
    return _json(find_cohorts(store(), Criteria.from_dict(criteria), limit=limit))


@app.get("/api/project/{accession}")
def project(accession: str):
    res = project_summary(store(), accession.upper())
    if "error" in res:
        raise HTTPException(404, res["error"])
    return _json(res)


@app.get("/api/files/{accession}")
def files(accession: str, where: str | None = None):
    res = sample_file_map(store(), accession.upper(), where)
    if "error" in res:
        raise HTTPException(404, res["error"])
    return _json(res)


@app.get("/api/manifest")
def manifest(acc: str = Query(..., description="comma-separated accessions"), q: str = "", format: str = "tsv"):
    man = build_manifest(store(), [a.strip().upper() for a in acc.split(",") if a.strip()], q)
    name = "cohort_manifest." + ("json" if format == "json" else "tsv")
    body = to_json(man) if format == "json" else to_tsv(man)
    return PlainTextResponse(body, headers={"Content-Disposition": f'attachment; filename="{name}"'})


@app.get("/api/map")
def project_map():
    path = INDEX_DIR / "map.json"
    if not path.exists():
        raise HTTPException(404, "project map not built; run `msagent map`")
    return FileResponse(path, media_type="application/json")


@app.get("/api/benchmark")
def benchmark():
    path = ROOT / "benchmarks" / "results.json"
    if not path.exists():
        raise HTTPException(404, "no benchmark results; run `msagent bench`")
    return FileResponse(path, media_type="application/json")

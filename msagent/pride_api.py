"""Optional enrichment from the PRIDE Archive REST API (v3).

Adds what SDRF files do not carry: project title, description, protocols, keywords,
submission/publication dates, publications, and the archive's own file listing (with
categories, sizes and download URLs). Responses are cached on disk.

Endpoints (same ones pridepy uses):
  GET {api}/projects/{accession}
  GET {api}/projects/{accession}/files/all
  GET {api}/search/projects?keyword=...&pageSize=...&page=...
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import duckdb
import httpx

from .config import CACHE_DIR, DB_PATH, PRIDE_API


class PrideUnavailable(RuntimeError):
    pass


def _client() -> httpx.Client:
    return httpx.Client(timeout=60, headers={"Accept": "application/json"}, follow_redirects=True)


def _get_json(client: httpx.Client, url: str, retries: int = 3, params: dict | None = None):
    last = None
    for attempt in range(retries):
        try:
            r = client.get(url, params=params)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except (httpx.ConnectError, httpx.ProxyError) as exc:
            raise PrideUnavailable(f"Cannot reach PRIDE API at {url}: {exc}") from exc
        except httpx.HTTPError as exc:
            last = exc
            time.sleep(2**attempt)
    raise PrideUnavailable(f"PRIDE API request failed for {url}: {last}")


def check_reachable() -> bool:
    try:
        with _client() as c:
            _get_json(c, f"{PRIDE_API}/projects/PXD000001", retries=1)
        return True
    except PrideUnavailable:
        return False


def fetch_project(accession: str, client: httpx.Client | None = None, cache_dir: Path = CACHE_DIR) -> dict | None:
    path = Path(cache_dir) / "pride" / f"{accession}.json"
    if path.exists():
        return json.loads(path.read_text())
    own = client is None
    client = client or _client()
    try:
        project = _get_json(client, f"{PRIDE_API}/projects/{accession}")
        files = _get_json(client, f"{PRIDE_API}/projects/{accession}/files/all") if project else None
    finally:
        if own:
            client.close()
    if project is None:
        return None
    record = {"project": project, "files": files or []}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record))
    return record


def search_projects(keyword: str, page_size: int = 100) -> list[str]:
    """PRIDE's own keyword search; returns accessions in PRIDE's ranking order."""
    with _client() as c:
        data = _get_json(
            c,
            f"{PRIDE_API}/search/projects",
            params={"keyword": keyword, "pageSize": page_size, "page": 0},
        )
    items = data
    if isinstance(data, dict):
        items = data.get("_embedded", {}).get("compactprojects") or data.get("content") or data.get("list") or []
    return [i.get("accession") for i in items or [] if isinstance(i, dict) and i.get("accession")]


def _names(items) -> list[str]:
    out = []
    for it in items or []:
        if isinstance(it, str):
            out.append(it)
        elif isinstance(it, dict):
            out.append(str(it.get("name") or it.get("value") or ""))
    return [o for o in out if o]


def parse_project(record: dict) -> tuple[tuple, list[tuple]]:
    p = record["project"]
    acc = p.get("accession")
    refs = p.get("references") or []
    pubmed = [str(r.get("pubmedID") or r.get("pubmedId")) for r in refs if r.get("pubmedID") or r.get("pubmedId")]
    dois = [str(r.get("doi")) for r in refs if r.get("doi")]
    kw = p.get("keywords") or []
    if isinstance(kw, str):
        kw = [k.strip() for k in kw.split(",") if k.strip()]
    prow = (
        acc,
        p.get("title"),
        p.get("projectDescription") or p.get("description"),
        p.get("sampleProcessingProtocol"),
        p.get("dataProcessingProtocol"),
        _names(kw) if kw and not isinstance(kw[0], str) else kw,
        (p.get("submissionDate") or "")[:10] or None,
        (p.get("publicationDate") or "")[:10] or None,
        pubmed,
        dois,
        json.dumps(p)[:200000],
    )
    frows = []
    for f in record.get("files") or []:
        locs = f.get("publicFileLocations") or []
        url = None
        for loc in locs:
            v = loc.get("value") or ""
            if v.startswith("ftp://ftp.pride.ebi.ac.uk/"):
                v = "https://ftp.pride.ebi.ac.uk/" + v[len("ftp://ftp.pride.ebi.ac.uk/"):]
            if v.startswith("http") or v.startswith("ftp"):
                url = v
                break
        cat = f.get("fileCategory")
        frows.append(
            (
                acc,
                f.get("fileName"),
                cat.get("value") if isinstance(cat, dict) else cat,
                f.get("fileSizeBytes"),
                url,
                f.get("checksum"),
            )
        )
    return prow, frows


def enrich(accessions: list[str] | None = None, db_path: Path = DB_PATH, limit: int | None = None, verbose: bool = True) -> dict:
    """Fetch PRIDE metadata for accessions in the DB and store it. Safe to re-run (cached)."""
    con = duckdb.connect(str(db_path))
    if accessions is None:
        accessions = [r[0] for r in con.execute(
            "SELECT accession FROM projects WHERE accession LIKE 'PXD%' AND accession NOT IN (SELECT accession FROM pride_projects) ORDER BY accession"
        ).fetchall()]
    if limit:
        accessions = accessions[:limit]
    if accessions and not check_reachable():
        con.close()
        raise PrideUnavailable(
            f"{PRIDE_API} is not reachable from this machine. Allow www.ebi.ac.uk in your network policy and re-run."
        )
    done = missing = 0
    with _client() as client:
        for i, acc in enumerate(accessions):
            rec = fetch_project(acc, client)
            if rec is None:
                missing += 1
                continue
            prow, frows = parse_project(rec)
            con.execute("DELETE FROM pride_projects WHERE accession = ?", [acc])
            con.execute("DELETE FROM pride_files WHERE accession = ?", [acc])
            con.execute("INSERT INTO pride_projects VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", list(prow))
            if frows:
                con.executemany("INSERT INTO pride_files VALUES (?, ?, ?, ?, ?, ?)", frows)
            done += 1
            if verbose and (i + 1) % 100 == 0:
                print(f"  enriched {i + 1}/{len(accessions)}", flush=True)
    con.close()
    return {"enriched": done, "not_found": missing}

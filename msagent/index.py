"""Text indexes over project documents: BM25 (keyword) and dense vectors (semantic)."""

from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path

import bm25s
import duckdb
import numpy as np

from .config import DB_PATH, EMBEDDER, INDEX_DIR
from .embed import get_embedder


def project_documents(con: duckdb.DuckDBPyConnection) -> tuple[list[str], list[str]]:
    """Per-project text: SDRF-derived summary plus PRIDE title/description/protocols when enriched."""
    rows = con.execute(
        """
        SELECT p.accession,
               concat_ws(' ', pp.title, pp.description, array_to_string(pp.keywords, ', '),
                         pp.sample_protocol, p.doc_text)
        FROM projects p LEFT JOIN pride_projects pp USING (accession)
        ORDER BY p.accession
        """
    ).fetchall()
    return [r[0] for r in rows], [r[1] or "" for r in rows]


def build(db_path: Path = DB_PATH, index_dir: Path = INDEX_DIR, embedder: str = EMBEDDER, verbose: bool = True) -> dict:
    t0 = time.time()
    con = duckdb.connect(str(db_path), read_only=True)
    accessions, docs = project_documents(con)
    con.close()
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)

    corpus_tokens = bm25s.tokenize(docs, stopwords="en", show_progress=False)
    retriever = bm25s.BM25()
    retriever.index(corpus_tokens, show_progress=False)
    retriever.save(str(index_dir / "bm25"))

    emb = get_embedder(embedder)
    vecs = []
    for i in range(0, len(docs), 2048):
        vecs.append(emb.embed([d[:4000] for d in docs[i:i + 2048]]))
    mat = np.vstack(vecs).astype(np.float32)
    np.save(index_dir / "embeddings.npy", mat)
    (index_dir / "accessions.json").write_text(json.dumps(accessions))
    meta = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "embedder": emb.name,
        "embedder_spec": embedder,
        "n_docs": len(docs),
        "dim": int(mat.shape[1]),
    }
    (index_dir / "meta.json").write_text(json.dumps(meta, indent=1))
    if verbose:
        print(f"  indexed {len(docs)} projects with BM25 + {emb.name} in {time.time() - t0:.1f}s")
    return meta


class TextIndex:
    def __init__(self, index_dir: Path = INDEX_DIR):
        index_dir = Path(index_dir)
        if not (index_dir / "meta.json").exists():
            raise FileNotFoundError(f"No search index at {index_dir}; run `msagent index`.")
        self.meta = json.loads((index_dir / "meta.json").read_text())
        self.accessions: list[str] = json.loads((index_dir / "accessions.json").read_text())
        self.pos = {a: i for i, a in enumerate(self.accessions)}
        self.bm25 = bm25s.BM25.load(str(index_dir / "bm25"))
        self.vectors = np.load(index_dir / "embeddings.npy")
        self._embedder = None

    @property
    def embedder(self):
        if self._embedder is None:
            self._embedder = get_embedder(self.meta.get("embedder_spec", EMBEDDER))
        return self._embedder

    def keyword(self, query: str, k: int = 100) -> list[tuple[str, float]]:
        if not query.strip():
            return []
        toks = bm25s.tokenize([query], stopwords="en", show_progress=False)
        if not toks.vocab:
            return []
        k = min(k, len(self.accessions))
        docs, scores = self.bm25.retrieve(toks, k=k, show_progress=False)
        return [(self.accessions[int(d)], float(s)) for d, s in zip(docs[0], scores[0]) if s > 0]

    def semantic(self, query: str, k: int = 100) -> list[tuple[str, float]]:
        if not query.strip():
            return []
        q = self.embedder.embed([query])[0]
        sims = self.vectors @ q
        top = np.argsort(-sims)[:k]
        return [(self.accessions[int(i)], float(sims[i])) for i in top]

    def hybrid(self, query: str, k: int = 100, rrf_k: int = 60) -> list[tuple[str, float]]:
        """Reciprocal-rank fusion of keyword and semantic rankings."""
        scores: dict[str, float] = {}
        for ranking in (self.keyword(query, k), self.semantic(query, k)):
            for rank, (acc, _) in enumerate(ranking):
                scores[acc] = scores.get(acc, 0.0) + 1.0 / (rrf_k + rank + 1)
        return sorted(scores.items(), key=lambda x: -x[1])[:k]

"""Text embedders. Default is WordLlama, whose weights ship in the PyPI wheel (fully offline)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import EMBEDDER


class Embedder:
    name = "base"

    def embed(self, texts: list[str]) -> np.ndarray:  # pragma: no cover - interface
        raise NotImplementedError


class WordLlamaEmbedder(Embedder):
    def __init__(self):
        import wordllama
        from wordllama import WordLlama

        # The wheel bundles weights + tokenizer; point the loader at the package dir so nothing is downloaded.
        self.model = WordLlama.load(cache_dir=Path(wordllama.__file__).parent, disable_download=True)
        self.name = "wordllama/l2_supercat_256"

    def embed(self, texts: list[str]) -> np.ndarray:
        v = np.asarray(self.model.embed(texts, norm=True), dtype=np.float32)
        return v


class FastEmbedEmbedder(Embedder):
    def __init__(self, model: str):
        from fastembed import TextEmbedding

        self.model = TextEmbedding(model)
        self.name = f"fastembed/{model}"

    def embed(self, texts: list[str]) -> np.ndarray:
        v = np.asarray(list(self.model.embed(texts)), dtype=np.float32)
        return v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-9)


_CACHE: dict[str, Embedder] = {}


def get_embedder(spec: str = EMBEDDER) -> Embedder:
    if spec not in _CACHE:
        if spec.startswith("fastembed:"):
            _CACHE[spec] = FastEmbedEmbedder(spec.split(":", 1)[1])
        else:
            _CACHE[spec] = WordLlamaEmbedder()
    return _CACHE[spec]

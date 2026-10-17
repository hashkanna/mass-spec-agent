"""Paths and settings. Everything can be overridden with environment variables."""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = Path(os.environ.get("MSAGENT_DATA_DIR", ROOT / "data"))

# Community SDRF annotations (https://github.com/bigbio/sdrf-annotated-datasets).
SDRF_REPO_URL = "https://github.com/bigbio/sdrf-annotated-datasets"
SDRF_DIR = Path(os.environ.get("MSAGENT_SDRF_DIR", DATA_DIR / "raw" / "sdrf-annotated-datasets"))

DB_PATH = Path(os.environ.get("MSAGENT_DB", DATA_DIR / "msagent.duckdb"))
INDEX_DIR = Path(os.environ.get("MSAGENT_INDEX_DIR", DATA_DIR / "index"))
CACHE_DIR = Path(os.environ.get("MSAGENT_CACHE_DIR", DATA_DIR / "cache"))

PRIDE_API = os.environ.get("MSAGENT_PRIDE_API", "https://www.ebi.ac.uk/pride/ws/archive/v3")
PRIDE_PROJECT_PAGE = "https://www.ebi.ac.uk/pride/archive/projects/{accession}"

# "wordllama" works fully offline (weights ship inside the PyPI wheel).
# "fastembed:<model>" gives better semantic recall when Hugging Face is reachable.
EMBEDDER = os.environ.get("MSAGENT_EMBEDDER", "wordllama")

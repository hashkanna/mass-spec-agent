#!/usr/bin/env bash
# Fetch the community-curated SDRF corpus (sample -> data file annotations for ~12k ProteomeXchange datasets).
# Source: https://github.com/bigbio/sdrf-annotated-datasets (Apache-2.0). ~750 MB checked out.
set -euo pipefail
DEST="${MSAGENT_SDRF_DIR:-$(dirname "$0")/../data/raw/sdrf-annotated-datasets}"
if [ -d "$DEST/.git" ]; then
  git -C "$DEST" pull --ff-only
else
  mkdir -p "$(dirname "$DEST")"
  git clone --depth 1 https://github.com/bigbio/sdrf-annotated-datasets "$DEST"
fi
git -C "$DEST" log -1 --format='SDRF corpus at %H (%cI)'

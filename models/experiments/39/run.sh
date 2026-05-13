#!/bin/bash
# Exp 39: disable noise filter (cnv_max_transitions=null) using exp 38 artifacts.
# No retrain, no HMM re-run — train.py skips both when reconstructions.npy and segments.parquet exist.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP38_DIR="$REPO_ROOT/data/results/38_06convvae_400bp_sinloss2.0_warmup80_eval05"
OUT_DIR="$REPO_ROOT/data/results/39_06convvae_400bp_sinloss2.0_warmup80_notransfilter"

mkdir -p "$OUT_DIR"

# Symlink large artifacts from exp 38 — skip training, inference, and HMM segmentation.
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy segments.parquet; do
    ln -sf "$EXP38_DIR/$f" "$OUT_DIR/$f"
done

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/39/config.yaml

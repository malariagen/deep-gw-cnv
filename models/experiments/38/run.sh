#!/bin/bash
# Exp 38: upgrade to 05_pf9_evaluation using exp 37 artifacts unchanged.
# No retrain, no HMM re-run — train.py skips both when reconstructions.npy and segments.parquet exist.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP37_DIR="$REPO_ROOT/data/results/37_06convvae_400bp_sinloss2.0_warmup80"
OUT_DIR="$REPO_ROOT/data/results/38_06convvae_400bp_sinloss2.0_warmup80_eval05"

mkdir -p "$OUT_DIR"

# Symlink large artifacts from exp 37 — skip training, inference, and HMM segmentation.
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy segments.parquet; do
    ln -sf "$EXP37_DIR/$f" "$OUT_DIR/$f"
done

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/38/config.yaml

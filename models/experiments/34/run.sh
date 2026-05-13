#!/bin/bash
# Exp 34: noise filter (v16) + CRR rescue on exp 32 artifacts. No retraining.
# Symlinks checkpoint, reconstructions, and segments from exp 32 so train.py
# skips straight to CNV calling with the noise filter and rescue threshold.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP32_DIR="$REPO_ROOT/data/results/32_06convvae_400bp_core_st0.90_hmm_primary_conf0.25"
OUT_DIR="$REPO_ROOT/data/results/34_06convvae_400bp_core_st0.90_noise_filter_crr_rescue"

mkdir -p "$OUT_DIR"

# Symlink all large artifacts from exp 32 — no retrain, no new HMM run needed.
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy segments.parquet; do
    ln -sf "$EXP32_DIR/$f" "$OUT_DIR/$f"
done

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/34/config.yaml

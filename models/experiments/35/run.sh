#!/bin/bash
# Exp 35: lower CRR rescue threshold (1.60 → 1.45) targeting GCH1 FNs in exp 34.
# No retraining, no new HMM run. Symlinks all large artifacts from exp 32.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP32_DIR="$REPO_ROOT/data/results/32_06convvae_400bp_core_st0.90_hmm_primary_conf0.25"
OUT_DIR="$REPO_ROOT/data/results/35_06convvae_400bp_core_st0.90_noise_filter_crr1.45"

mkdir -p "$OUT_DIR"

# Symlink all large artifacts from exp 32 — no retrain, no new HMM run needed.
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy segments.parquet; do
    ln -sf "$EXP32_DIR/$f" "$OUT_DIR/$f"
done

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/35/config.yaml

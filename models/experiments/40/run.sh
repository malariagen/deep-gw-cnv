#!/bin/bash
# Exp 40: raise hmm_self_transition 0.90 → 0.96 to fix sanity-check mass failure.
# Reuses exp 37 checkpoint, latents, reconstructions. Re-runs HMM (segments.parquet
# NOT symlinked), then CNV calling and evaluation.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP37_DIR="$REPO_ROOT/data/results/37_06convvae_400bp_sinloss2.0_warmup80"
OUT_DIR="$REPO_ROOT/data/results/40_06convvae_400bp_sinloss2.0_warmup80_st096"

mkdir -p "$OUT_DIR"

# Symlink VAE artifacts — skip training and inference. Do NOT symlink segments.parquet;
# HMM must re-run with the new self_transition=0.96.
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy; do
    ln -sf "$EXP37_DIR/$f" "$OUT_DIR/$f"
done

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/40/config.yaml

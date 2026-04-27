#!/bin/bash
# Exp 22: VAE fine-tune from exp 21 checkpoint with synthetic normal augmentation.
# aug_normal_poisson=true Poisson-resamples normal samples each epoch to tighten
# the VAE's normal prior, aiming to raise FN CRR above the 1.30 gate without
# lowering the gate. Full retraining required (no skip-training path).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/22/config.yaml

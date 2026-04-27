#!/bin/bash
# Exp 21: curriculum-weighted VAE fine-tune from exp 18 checkpoint (final ratio 0.30).
# Runs full training (VAE + HMM + CNV calling + evaluation).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/21/config.yaml

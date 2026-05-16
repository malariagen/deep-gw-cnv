#!/bin/bash
# Exp 42: retrain 06_conv_vae without sin loss. Full train + HMM + CNV calling + evaluation.
# Sin loss (max_weight=2.0) produced noisy reconstructions; removing it requires a retrain.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/42/config.yaml

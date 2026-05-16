#!/bin/bash
# Exp 43: retrain 06_conv_vae with max_beta=4.0 + 18_genome_cnv_caller (rescue before sanity).
# Full train + HMM + CNV calling + evaluation.
# max_beta change alters the training objective, so exp 42 checkpoint cannot be reused.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/43/config.yaml

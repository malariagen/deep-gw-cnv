#!/bin/bash
# Exp 45: retrain 06_conv_vae with beta=4.0, cnv_downsample_ratio_final=0.05 + 18_genome_cnv_caller.
# Full train + HMM + CNV calling + evaluation.
# Curriculum change (ramp eliminated; flat 5% throughout) requires full retrain.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/45/config.yaml

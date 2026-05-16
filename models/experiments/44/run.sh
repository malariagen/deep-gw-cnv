#!/bin/bash
# Exp 44: retrain 06_conv_vae with beta=4.0, cnv_downsample_ratio_final=0.10 + 18_genome_cnv_caller.
# Full train + HMM + CNV calling + evaluation.
# Curriculum change alters training dynamics; exp 43 checkpoint cannot be reused.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/44/config.yaml

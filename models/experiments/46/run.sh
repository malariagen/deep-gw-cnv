#!/bin/bash
# Exp 46 (revised): fully-connected VAE (07_fc_vae), beta=4.0, flat 5% curriculum, 18_genome_cnv_caller.
# Training completed in the original run. This script resumes from the checkpoint via wrap_up.py,
# skipping training and running inference → HMM → CNV calling → evaluation only.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" -m training.wrap_up experiments/46/config.yaml

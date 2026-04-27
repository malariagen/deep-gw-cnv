#!/bin/bash
# Exp 29: 06_conv_vae on 400 bp core genome — full training from scratch.
# No checkpoint reuse (linear layer dimensions change: flat_dim 166656 → 416000).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/29/config.yaml

#!/bin/bash
# Exp 33: Deeper VAE training (200 epochs). Full training from scratch.
# Better reconstruction quality should improve HMM signal at gene loci.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

"$REPO_ROOT/.venv/bin/python" train.py experiments/33/config.yaml

#!/bin/bash
# Exp 37: full retrain of 06_conv_vae at 400bp with boosted sinusoidal loss and slower curriculum.
# Addresses GCH1 frac=0 failure in AF-W/AF-E by improving VAE reconstruction sensitivity.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" train.py experiments/37/config.yaml

#!/bin/bash
# Exp 32: lower HMM self_transition 0.96 → 0.90 to recover gene-scale FNs.
# Reuses exp 30 checkpoint, latents, reconstructions, and sample_ids via symlinks.
# HMM must be re-run (new st value) — segments.parquet is NOT symlinked.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

OUT="$REPO_ROOT/data/results/32_06convvae_400bp_core_st0.90_hmm_primary_conf0.25"
SRC="$REPO_ROOT/data/results/30_06convvae_400bp_core_st0.96_cnvcaller_fix"

mkdir -p "$OUT"
# Reuse inference outputs (same checkpoint, same data)
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy; do
    ln -sf "$SRC/$f" "$OUT/$f"
done

"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, yaml, importlib, pathlib
sys.path.insert(0, ".")

cfg_path = "experiments/32/config.yaml"
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)

resolve = lambda p: str((pathlib.Path(cfg_path).parent / p).resolve())
store_path = resolve(cfg["store_path"])
out_dir    = resolve(cfg["out_dir"])

run_hmm_all_samples = importlib.import_module(f"hmm.{cfg['hmm']}").run_hmm_all_samples
run_cnv_calls       = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
run_evaluation      = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

cfg_r = dict(cfg)
cfg_r["gff_path"]      = resolve(cfg["gff_path"])
cfg_r["pf9_gt_path"]   = resolve(cfg["pf9_gt_path"])
cfg_r["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])

print("Running HMM segmentation (st=0.90)...", flush=True)
run_hmm_all_samples(store_path, out_dir, cfg_r)

print("Calling genome-wide CNVs...", flush=True)
run_cnv_calls(store_path, out_dir, cfg_r)

print("Running evaluation...", flush=True)
run_evaluation(out_dir, cfg_r)
print("Done.")
EOF

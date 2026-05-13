#!/bin/bash
# Exp 31: HMM-segment-primary calling (v13 caller; CRR diagnostic only).
# Reuses exp 30 checkpoint, latents, reconstructions, sample_ids, and segments via symlinks.
# Only CNV calling and evaluation re-run — no HMM needed.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

OUT="$REPO_ROOT/data/results/31_06convvae_400bp_core_st0.96_hmm_primary_conf0.25"
SRC="$REPO_ROOT/data/results/30_06convvae_400bp_core_st0.96_cnvcaller_fix"

mkdir -p "$OUT"
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy segments.parquet; do
    ln -sf "$SRC/$f" "$OUT/$f"
done

"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, yaml, importlib, pathlib
sys.path.insert(0, ".")

cfg_path = "experiments/31/config.yaml"
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)

resolve = lambda p: str((pathlib.Path(cfg_path).parent / p).resolve())
store_path = resolve(cfg["store_path"])
out_dir    = resolve(cfg["out_dir"])

run_cnv_calls  = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
run_evaluation = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

cfg_r = dict(cfg)
cfg_r["gff_path"]      = resolve(cfg["gff_path"])
cfg_r["pf9_gt_path"]   = resolve(cfg["pf9_gt_path"])
cfg_r["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])

print("Calling genome-wide CNVs...", flush=True)
run_cnv_calls(store_path, out_dir, cfg_r)

print("Running evaluation...", flush=True)
run_evaluation(out_dir, cfg_r)
print("Done.")
EOF

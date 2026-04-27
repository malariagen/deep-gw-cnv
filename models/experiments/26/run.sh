#!/bin/bash
# Exp 26: Band segment filter — demote CRR [1.20, 1.25) calls not backed by CN>=2 HMM segments.
# Reuses exp 25 checkpoint, latents, reconstructions, sample_ids, AND segments via symlinks.
# HMM is unchanged (self_transition=0.90 from exp 25); only CNV caller and evaluation re-run.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

OUT="$REPO_ROOT/data/results/26_skip_st0.90_sin_curriculum_r0.30_genome_wide_band_segcov0.5"

mkdir -p "$OUT"
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy segments.parquet; do
    ln -sf "../25_skip_st0.90_sin_curriculum_r0.30_genome_wide_thr1.20/$f" "$OUT/$f"
done

"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, yaml, importlib, pathlib
sys.path.insert(0, ".")

cfg_path = "experiments/26/config.yaml"
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)

resolve = lambda p: str((pathlib.Path(cfg_path).parent / p).resolve())
store_path = resolve(cfg["store_path"])
out_dir    = resolve(cfg["out_dir"])

run_cnv_calls  = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
run_evaluation = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

cfg_resolved = dict(cfg)
cfg_resolved["gff_path"]      = resolve(cfg["gff_path"])
cfg_resolved["pf9_gt_path"]   = resolve(cfg["pf9_gt_path"])
cfg_resolved["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])

print("Calling genome-wide CNVs (band segment filter: CRR [1.20, 1.25) needs CN>=2 coverage >= 0.5)...", flush=True)
run_cnv_calls(store_path, out_dir, cfg_resolved)

print("Running evaluation on reference genes...", flush=True)
run_evaluation(out_dir, cfg_resolved)
print("Done.")
EOF

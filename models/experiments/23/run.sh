#!/bin/bash
# Exp 23: Genome-wide CNV calling across all P. falciparum genes, reusing exp 21 artifacts.
# No retraining or HMM re-run: latents, reconstructions, sample_ids, segments,
# and checkpoint are symlinked from exp 21. Only CNV calling and eval are re-run.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

EXP21="$REPO_ROOT/data/results/21_skip_st0.80_sin_curriculum_r0.30"
OUT="$REPO_ROOT/data/results/23_skip_st0.80_sin_curriculum_r0.30_genome_wide"

mkdir -p "$OUT"
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy segments.parquet; do
    ln -sf "../21_skip_st0.80_sin_curriculum_r0.30/$f" "$OUT/$f"
done

"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, yaml, importlib, pathlib
sys.path.insert(0, ".")

cfg_path = "experiments/23/config.yaml"
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

print("Calling genome-wide CNVs...", flush=True)
run_cnv_calls(store_path, out_dir, cfg_resolved)

print("Running evaluation on reference genes...", flush=True)
run_evaluation(out_dir, cfg_resolved)
print("Done.")
EOF

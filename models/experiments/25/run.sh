#!/bin/bash
# Exp 25: Raise HMM self_transition 0.80→0.90 to reduce jumpy micro-calls and FPs.
# Reuses exp 21 checkpoint, latents, reconstructions, sample_ids via symlinks.
# HMM is re-run (not symlinked) because self_transition changed.
# CNV calling and evaluation re-run with same 08 caller and 1.20 threshold.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

OUT="$REPO_ROOT/data/results/25_skip_st0.90_sin_curriculum_r0.30_genome_wide_thr1.20"

mkdir -p "$OUT"
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy; do
    ln -sf "../21_skip_st0.80_sin_curriculum_r0.30/$f" "$OUT/$f"
done

"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, yaml, importlib, pathlib
sys.path.insert(0, ".")

cfg_path = "experiments/25/config.yaml"
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)

resolve = lambda p: str((pathlib.Path(cfg_path).parent / p).resolve())
store_path = resolve(cfg["store_path"])
out_dir    = resolve(cfg["out_dir"])

run_hmm_all_samples = importlib.import_module(f"hmm.{cfg['hmm']}").run_hmm_all_samples
run_cnv_calls       = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
run_evaluation      = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

cfg_resolved = dict(cfg)
cfg_resolved["gff_path"]      = resolve(cfg["gff_path"])
cfg_resolved["pf9_gt_path"]   = resolve(cfg["pf9_gt_path"])
cfg_resolved["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])

print("Fitting HMM segments (self_transition=0.90)...", flush=True)
run_hmm_all_samples(store_path, out_dir, cfg)

print("Calling genome-wide CNVs...", flush=True)
run_cnv_calls(store_path, out_dir, cfg_resolved)

print("Running evaluation on reference genes...", flush=True)
run_evaluation(out_dir, cfg_resolved)
print("Done.")
EOF

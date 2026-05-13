#!/bin/bash
# Exp 30: fix CNV caller CRR bug + recalibrated HMM self_transition (0.96) for 400bp.
# Reuses exp 29 checkpoint, latents, reconstructions, sample_ids via symlinks.
# Only HMM, CNV calling, and evaluation are re-run.
#
# HMM and CNV caller run in separate Python processes so the HMM's working memory
# (joblib workers, mmap handles) is fully released before the CNV caller loads the
# full counts and reconstructions arrays.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

OUT="$REPO_ROOT/data/results/30_06convvae_400bp_core_st0.96_cnvcaller_fix"
SRC="$REPO_ROOT/data/results/29_06convvae_400bp_core_st0.90_sin_curriculum_r0.30"

mkdir -p "$OUT"
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy; do
    ln -sf "$SRC/$f" "$OUT/$f"
done

# ── Step 1: HMM segmentation ────────────────────────────────────────────────
"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, yaml, importlib, pathlib
sys.path.insert(0, ".")

cfg_path = "experiments/30/config.yaml"
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)

resolve = lambda p: str((pathlib.Path(cfg_path).parent / p).resolve())

run_hmm_all_samples = importlib.import_module(f"hmm.{cfg['hmm']}").run_hmm_all_samples

cfg_r = dict(cfg)
print("Fitting HMM (self_transition=0.96)...", flush=True)
run_hmm_all_samples(resolve(cfg["store_path"]), resolve(cfg["out_dir"]), cfg_r)
print("HMM done.", flush=True)
EOF

# ── Step 2: CNV calling + evaluation ────────────────────────────────────────
"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, yaml, importlib, pathlib
sys.path.insert(0, ".")

cfg_path = "experiments/30/config.yaml"
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

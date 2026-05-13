# Experiment 30 — CNV caller CRR bug fix + HMM self_transition recalibration

**Status:** Complete 2026-04-28

**No retraining.** Reuses exp 29 checkpoint, latents, and reconstructions via symlinks.
Only HMM segmentation, CNV calling, and evaluation are re-run.

## Hypothesis

Exp 29 had two independent problems.

**Problem 1 — CNV caller arithmetic bug (primary; explains MDR1 FNR=0.88).**
In `11_genome_cnv_caller.py`, CRR is computed as `mean(gene copy_ratio) / mean(flank copy_ratio)`.
The copy_ratio per bin is `counts / (recons + 1e-6)`. At 400bp resolution, ~2% of bins
per sample have reconstruction near zero while carrying non-zero read counts. Each such
bin produces a copy_ratio in the millions. With even 2–5 such bins in the flank (out of
~2689 on chr5), `nanmean(flank copy_ratio)` collapses to thousands, making CRR ≈ 0 for
almost every sample regardless of true amplification status.

This is a numerical artifact, not a model quality issue. Visual inspection in the
Streamlit app confirmed the reconstruction does track the input at MDR1. Retrospective
re-computation with low-coverage bins excluded (counts < 10 OR recons < 10) gives
MDR1 GT-amplified CRR p50 = 2.14 (well above the 1.20 threshold) and GT-normal
CRR p50 = 1.01 — perfectly separated.

**Problem 2 — HMM over-segmentation (secondary; explains CRT and PM2_PM3 regressions).**
At 400bp, self_transition=0.90 produces p50=270 segment transitions per sample vs p50=62
at 1000bp, fragmenting gene-spanning CN=2 runs and causing the minimum-gene-coverage
filter to block valid calls. Scaling to target_st = 1 − (1−0.90)×(1500/3700) ≈ 0.96
restores equivalent coarseness.

## Changes from experiment 29

| Component              | Exp 29 | Exp 30 | Rationale                                                         |
|------------------------|--------|--------|-------------------------------------------------------------------|
| `cnv`                  | `11_genome_cnv_caller` | **`12_genome_cnv_caller`** | Fixes near-zero-recon flank contamination |
| `hmm_self_transition`  | 0.90   | **0.96** | Calibrated for 2.5× more bins/chromosome at 400bp              |
| Checkpoint             | trained fresh | reused from exp 29 | No model change requires retraining            |
| Inference outputs      | computed | symlinked from exp 29 | Identical model → identical outputs          |

`12_genome_cnv_caller` is identical to v11 except it NaN-masks bins where
`counts < hmm_low_cov_threshold OR recons < hmm_low_cov_threshold` before computing
the gene and flank means. `np.nanmean` then excludes those bins automatically.

## Expected outcome

- **MDR1 FNR**: 0.88 → ~0.05–0.20. Corrected CRR has strong signal (GT-amplified
  p50=2.14 before HMM gate); remaining FNs will depend on HMM coverage of the MDR1
  segment, which should now be less fragmented with st=0.96.
- **CRT FNR**: 0.22 → ~0.05–0.15. HMM fix reduces fragmentation; CRR for CRT FNs
  (p50=0.00 in exp 29) was also contaminated by the same caller bug.
- **PM2_PM3 FNR**: 0.14 → ~0.05–0.10. Same dual-fix benefit.
- **GCH1 FNR**: 0.08 → ~0.08. GCH1 was not affected by either issue; no regression expected.

## Run time estimate

~2–3 hours: HMM (~2h with n_jobs=-1 across 51986 bins × 54k samples) +
CNV calling (~30 min) + evaluation (~5 min). No training.

## Actual outcome

| Gene    | MCC  | FNR  | PPV  |
|---------|------|------|------|
| CRT     | 0.27 | 0.02 | 0.07 |
| GCH1    | 0.75 | 0.03 | 0.63 |
| MDR1    | 0.71 | 0.00 | 0.52 |
| PM2_PM3 | 0.69 | 0.02 | 0.50 |

**Where predictions were right:** FNR improvements exceeded expectations across all
genes. MDR1 FNR dropped to 0.00 (predicted 0.05–0.20); CRT to 0.02 (predicted
0.05–0.15); PM2_PM3 to 0.02 (predicted 0.05–0.10); GCH1 to 0.03 (predicted ~0.08).
The caller bug fix and HMM recalibration together worked as hypothesised.

**Where predictions diverged — PPV collapsed:** The move to 400bp plus the caller
fix surfaced a large number of low-CRR positive calls. CRT PPV=0.07 (was 0.42 at
1000bp); MDR1 PPV=0.52 (was 0.91). CRR by outcome shows FPs concentrate at
CRR p50≈1.40, while TPs sit at p50≈2.0–2.3. The FPs with CRR 1.35–1.55 are
currently passing unconditionally because they exceed the `cnv_crr_band_upper=1.35`
threshold. These calls bypass the HMM CN=2 coverage filter entirely.

**Root cause of low PPV:** `cnv_crr_band_upper=1.35` (the unconditional-pass
boundary) is too low at 400bp. The more numerous bins per gene/flank give noisier
CRR estimates, so some normal samples land in the 1.35–1.55 range. At 1000bp
fewer bins meant less noise and fewer false exceedances of 1.35.

**Next experiment:** Exp 31 — raise `cnv_crr_band_upper` from 1.35 → 1.55,
gating the 1.35–1.55 range through the 20% CN=2 coverage check. No HMM re-run
needed; reuses exp 30 segments.

## Proposal history

**Original proposal (2026-04-28):** MDR1 FNR=0.88 was attributed to the VAE absorbing
MDR1 amplification as "normal" (CRR p50=0.00 for FNs interpreted as reconstruction >>
input at MDR1). The original proposal was to fix only the HMM: self_transition 0.90 →
0.96. MDR1 was left unsolved, with a separate curriculum-correction retraining flagged
as the intended follow-up.

**Feedback received (2026-04-28):** MDR1 inspected in the Streamlit app appeared fine —
the reconstruction tracks the input well at MDR1. Requested investigation of the
discrepancy between visual inspection and the evaluation numbers.

**What changed:** Investigation revealed the CRR collapse was caused by a bug in
`11_genome_cnv_caller.py`, not the model. Near-zero reconstruction bins in the flank
inflate the flank mean by factors of thousands. The fix (v12 caller) excludes these
bins. MDR1 is now expected to recover strongly without any retraining.

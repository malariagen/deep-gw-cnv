# Experiment 36 — Global CRR threshold 1.55 + population diagnostic

**Status:** Complete 2026-05-06

**No retraining, no new HMM run.** Symlinks exp 32 checkpoint, reconstructions, and
segments. Single config change from exp 35: `cnv_crr_rescue_threshold` 1.45 → 1.55.
Reverts from the originally-proposed per-gene v17 caller.

## Hypothesis

Exp 35 achieved GCH1 FNR=0.09 by lowering the global CRR threshold to 1.45, but this
introduced two regressions: CRT PPV 0.92 → 0.55 (16 new FPs with CRR 1.45–1.52) and
PM2_PM3 PPV 0.90 → 0.81 (~30 new FPs with CRR 1.46–1.63).

The CNV caller is genome-wide and must use a single global threshold — per-gene thresholds
would not generalise beyond the four reference genes. The compromise threshold is 1.55:
three basis points above the CRT FP ceiling (p90=1.52), which eliminates all 16 CRT FPs
while still rescuing GCH1 and MDR1 FNs with CRR ≥ 1.55.

The secondary purpose is to establish a stable baseline for population-level diagnostics.
The evaluation BY POPULATION breakdown will show whether the remaining GCH1 FNR is
concentrated in AF-W and AF-E — if so, the residual problem is a population-specific
signal issue, not a threshold issue, and the next experiment should target the model rather
than the calling parameters.

## Changes from experiment 35

| Parameter                  | Exp 35  | Exp 36   | Rationale                                                  |
|----------------------------|---------|----------|------------------------------------------------------------|
| `cnv_crr_rescue_threshold` | 1.45    | **1.55** | Eliminates all CRT FPs (ceiling 1.52); partial PM2_PM3 fix |

All other parameters and all artifacts unchanged.

## Expected outcome

- **CRT PPV**: 0.55 → ~0.92+. All 16 FPs have CRR ≤ 1.52; threshold 1.55 eliminates them
  without touching TPs (all CRR ≥ 1.66).
- **PM2_PM3 PPV**: 0.81 → ~0.84–0.87. FPs with CRR in [1.46, 1.54] are eliminated;
  those with CRR ≥ 1.55 remain.
- **GCH1 FNR**: 0.09 → ~0.16–0.18. The 66 remaining FNs (CRR ≤ 1.43) are unaffected.
  The 81 FNs rescued at threshold 1.45 with CRR in [1.45, 1.54] become FN again; only
  those with CRR ≥ 1.55 remain rescued (~15–20% of the original 147, ≈ 25–30 FNs).
- **MDR1**: FNR already 0.01; minimal effect.
- **Missingness and call_rate**: unchanged (noise filter unchanged).

## Population diagnostic

The evaluation BY POPULATION table will show GCH1 FNR per population group. The working
hypothesis is that AF-W and AF-E drive the residual FNR. If the BY POPULATION table
confirms high FNR concentrated in those two populations while AS-SE-E and OC-NG remain
near zero, the next experiment should investigate whether a different VAE architecture or
normalisation scheme produces better CRR signal in low-prevalence African populations.

## Run time estimate

~30 minutes (CNV calling + evaluation only; all artifacts symlinked from exp 32).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.04 | 0.92 | 0.94 |
| GCH1    | 0.16 | 0.98 | 0.90 |
| MDR1    | 0.02 | 0.89 | 0.93 |
| PM2_PM3 | 0.00 | 0.88 | 0.93 |

**Predictions matched:** CRT PPV recovered to 0.92 (all 16 FPs eliminated, as predicted). GCH1 FNR
regressed to 0.16 (within predicted 0.16–0.18 range). PM2_PM3 PPV reached 0.88, slightly better
than the predicted 0.84–0.87.

**Population diagnostic confirmed hypothesis:** GCH1 FNR is severely concentrated in AF-W (FNR=0.49)
and AF-E (FNR=0.20). AS-SE-E and OC-NG remain near 0. This is a population-specific signal failure,
not a threshold issue.

**Near-miss diagnostic:** 116/117 GCH1 FNs have segment_cn2_fraction=0 — the HMM produces no
elevated-CN coverage for those gene bins at all. FN CRR p50=1.42 confirms the raw signal is present
but the VAE reconstruction doesn't amplify it enough for HMM detection. Parameter tuning on calling
thresholds cannot fix frac=0 cases; the signal must be improved upstream (model or bin resolution).

**Next experiment:** Exp 37 — retrain 06_conv_vae at 400bp with boosted sinusoidal loss (1.0 → 2.0) and slower curriculum warmup (40 → 80 epochs).

## Proposal history

**First proposal (2026-05-06):** Experiment 36 was originally proposed as a new v17 caller
introducing a per-gene CRR rescue threshold dict (GCH1: 1.35, CRT: 1.55, MDR1: 1.45,
PM2_PM3: 1.50). This was designed to break the deadlock between GCH1 FNR (requiring a
low threshold) and CRT/PM2_PM3 PPV (requiring a high threshold).

**Feedback received:** Chiyun rejected the per-gene approach, noting that the goal is to
call CNVs genome-wide — a gene-keyed threshold dict only serves the four named reference
genes and does not generalise to any other gene. Chiyun approved the alternative direction
of population-level diagnostics to understand whether the residual GCH1 FNs are
structurally different in AF-W/AF-E populations.

**Revised proposal:** Drop v17 caller, revert to v16 with a single global threshold of
1.55. This fixes the CRT/PM2_PM3 regression at the cost of partial GCH1 FNR regression,
and positions the evaluation output as a population diagnostic.

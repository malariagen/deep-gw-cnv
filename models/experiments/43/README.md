# Experiment 43 — Retrain 06_conv_vae with β=4 + 18_genome_cnv_caller (rescue before sanity)

**Status:** Complete 2026-05-16

**Full retrain.** `max_beta` raised 1.0 → 4.0 changes the training objective; the exp 42
checkpoint cannot be reused. HMM segmentation and CNV calling are also re-run.

## Hypothesis

Experiment 42 (sin loss removed, β=1) reduced HMM transitions p50 from 769 → 306 and drove
evaluable FNR to zero for CRT, MDR1, and PM2_PM3. Two problems remain:

1. **Masked CNV pool:** total FNR is 0.12–0.33 because the sanity check in `16_genome_cnv_caller`
   runs before CRR rescue. Samples oscillating enough to fail sanity (< 55% CN=1) are assigned
   CN=-1 before the rescue criterion fires, even when their raw CRR >> 1.55. Fix: `18_genome_cnv_caller`
   elevates CRR rescue above the sanity check.

2. **GCH1 evaluable FNR=0.14 (244 FNs, CRR p50=1.43):** these samples are not masked —
   the HMM is calling them normal (frac=0) because the reconstruction residual is modest.
   A higher β strengthens KL regularisation, pushing the VAE toward a more structured latent
   space where normal profiles reconstruct well and amplified profiles reconstruct poorly,
   potentially raising CRR contrast for the 244 GCH1 FNs.

## Changes from experiment 42

| Parameter / Component       | Exp 42                  | Exp 43                      | Rationale                                                                    |
|-----------------------------|-------------------------|-----------------------------|------------------------------------------------------------------------------|
| `max_beta`                  | 1.0                     | **4.0**                     | Standard β-VAE value; higher KL weight may improve CRR contrast             |
| `cnv` module                | `16_genome_cnv_caller`  | **`18_genome_cnv_caller`**  | Rescue applied before sanity check; same global threshold                    |
| VAE retrain                 | yes                     | **yes** (new β)             | Training objective changed; full retrain required                            |
| HMM re-run                  | yes                     | **yes** (new checkpoint)    | New reconstructions require new segmentation                                 |

All other parameters unchanged from exp 42.

## Expected outcome

- **β=4 effect on normal samples:** tighter KL regularisation → lower reconstruction variance
  for normal read-count profiles → reduced HMM oscillation. Transitions p50 may drop further
  below 306.
- **β=4 effect on amplified samples:** the VAE will reconstruct elevated profiles less well
  → higher CRR → potentially rescuing the 244 GCH1 evaluable FNs above the 1.55 threshold.
- **Masked CNV pool:** `18_genome_cnv_caller` should rescue most of the 344 GCH1 / 135 MDR1
  / 310 PM2_PM3 masked positives (CRR >> 1.55). Total FNR should converge toward evaluable FNR.
- **Risk:** if β=4 causes posterior collapse, the reconstruction signal weakens for all samples.
  Monitor: callability and CRR distributions. If callability drops substantially below 0.995,
  β is too high.

## Run time estimate

~7–8 hours: VAE training (~5–6h for 100 epochs) + HMM re-segmentation (~2h) + CNV calling
(~30 min) + evaluation (~5 min).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.02 | 0.56 | 0.74 |
| GCH1    | 0.15 | 0.97 | 0.90 |
| MDR1    | 0.00 | 0.87 | 0.93 |
| PM2_PM3 | 0.02 | 0.89 | 0.93 |

Masked CNV: CRT=0, GCH1=1, MDR1=1, PM2_PM3=9 (masked rate ≈ 0% for all genes).
Total FNR: CRT=0.02, GCH1=0.15, MDR1=0.01, PM2_PM3=0.03 — now converged to evaluable FNR.
HMM transitions p50=127 (down from 306 in exp 42; further 57% reduction). Callability: 0.995.

**Where predictions matched:** The `18_genome_cnv_caller` fix eliminated the masked CNV pool
almost entirely (344/135/310 masked → 0/1/9), confirming the caller ordering was the blocker.
Total FNR converged to evaluable FNR as expected. Callability held at 0.995 — β=4 caused no
posterior collapse. HMM stability improved substantially (p50: 306→127).

**Where predictions diverged:** β=4 did not improve CRR contrast for the 306 GCH1 evaluable FNs
(CRR p50: 1.43→1.42, essentially unchanged). All 306 FNs have frac=0 (HMM gives no CN>=2
coverage) and CRR below the 1.55 rescue threshold — these are not rescuable by threshold
changes and cannot be fixed by HMM parameter tuning alone. Root cause: the VAE sees enough
GCH1-amplified profiles in training (30% CNV in steady-state, high African prevalence) that it
reconstructs amplified GCH1 as near-normal, keeping reconstruction residual low.

**Next:** Experiment 44 — lower `cnv_downsample_ratio_final` from 0.30 → 0.10 to reduce the
proportion of amplified profiles seen in late training, pushing the VAE to treat amplification
as genuinely anomalous.

## Proposal history

**Original proposal (2026-05-15):** Caller-only fix — `18_genome_cnv_caller` elevating CRR
rescue above the sanity check, reusing exp 42 checkpoint and segments via symlinks. No VAE
retrain. Rationale: masked CRR p10=1.58–1.62 confirmed signal was present; the bottleneck
was caller ordering, not reconstruction quality.

**Feedback:** Can we increase beta?

**What changed:** `max_beta` raised 1.0 → 4.0, requiring a full VAE retrain and HMM re-run.
The caller fix (`18_genome_cnv_caller`) is retained as it independently addresses the masked
CNV pool. The two changes combine: a better-regularised VAE may raise CRR contrast for the
244 GCH1 evaluable FNs, while the rescue-before-sanity ordering recovers the existing masked pool.

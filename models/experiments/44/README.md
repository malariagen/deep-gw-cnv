# Experiment 44 — Lower CNV training proportion: cnv_downsample_ratio_final 0.30 → 0.10

**Status:** Proposed 2026-05-16

**Full retrain.** Curriculum change alters training dynamics from epoch 80 onward; the exp 43
checkpoint cannot be reused. HMM segmentation and CNV calling are also re-run.

## Hypothesis

Experiment 43 (β=4, `18_genome_cnv_caller`) eliminated the masked CNV pool and drove
HMM transitions p50 from 306 → 127. One problem remains:

**GCH1 evaluable FNR=0.15 (306 FNs, CRR p50=1.42, all frac=0):** The HMM never segments
these samples as CN>=2. Their CRR is real (1.42 on average, well above baseline) but the
VAE reconstruction residual is too small to trigger the HMM. β=4 did not improve CRR contrast
for these FNs (p50 moved 1.43→1.42 vs exp 42).

The root cause is the curriculum-weighted sampler's steady-state CNV proportion of 30%.
GCH1 amplification is highly prevalent in African populations (often 40–60%+ in AF-E, AF-W,
AF-C), which are the dominant populations in Pf9. Even with downsampling to 30%, the VAE
sees enough amplified GCH1 profiles during training that it learns elevated read counts as
near-normal for that locus — keeping the reconstruction residual low and the HMM uninformed.

Lowering `cnv_downsample_ratio_final` to 10% makes amplification genuinely anomalous:
the VAE will see it rarely enough in late training that it reconstructs amplified profiles
poorly, raising residual and giving the HMM a signal to act on.

## Changes from experiment 43

| Parameter                        | Exp 43 | Exp 44    | Rationale                                                                              |
|----------------------------------|--------|-----------|----------------------------------------------------------------------------------------|
| `cnv_downsample_ratio_final`     | 0.30   | **0.10**  | Reduce CNV exposure in late training; force VAE to treat amplification as anomalous    |
| VAE retrain                      | yes    | **yes**   | Curriculum change alters training; checkpoint not reusable                             |
| HMM re-run                       | yes    | **yes**   | New reconstructions require new segmentation                                           |

All other parameters unchanged from exp 43 (β=4, `18_genome_cnv_caller`, global
`cnv_crr_rescue_threshold=1.55`).

## Expected outcome

- **GCH1 FNR:** With lower CNV exposure, the VAE should reconstruct amplified GCH1 profiles
  more poorly → higher reconstruction residual → HMM assigns CN>=2 → some frac=0 FNs become
  frac>0. If enough cross the 0.50 gene-coverage threshold, total GCH1 FNR could drop toward
  0.05–0.10.
- **HMM transitions p50:** May rise slightly (the VAE is less fitted to amplified profiles,
  so its residual variance may be higher for those samples). Monitor: if p50 climbs above 250,
  the trade-off may not be worthwhile.
- **Callability:** Should remain near 0.995; this change does not affect KL regularisation.
- **Other genes (CRT/MDR1/PM2_PM3):** FNR is already near-zero; small risk of PPV degradation
  if the VAE starts flagging more borderline profiles, but with only 10% CNV the model should
  still learn a clean normal-vs-amplified boundary.
- **Risk:** If 10% is too low, the VAE may become overly sensitive and produce noisy residuals
  for all samples, raising transitions p50. Key diagnostic: if transitions p50 > 250 and
  GCH1 FNR does not improve, the 10% ratio is too aggressive and the next step would be 0.15
  or 0.20.

## Run time estimate

~7–8 hours: VAE training (~5–6h for 100 epochs) + HMM re-segmentation (~2h) + CNV calling
(~30 min) + evaluation (~5 min).

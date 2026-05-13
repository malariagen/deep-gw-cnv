# Experiment 37 — Retrain 06_conv_vae at 400bp with boosted sinusoidal loss and slower curriculum

**Status:** Complete 2026-05-13

**Full retrain from scratch.** Two training-regimen changes to improve VAE reconstruction sensitivity
for weak African GCH1 amplification.

## Hypothesis

Experiment 36 confirmed that GCH1 FNR=0.16 is concentrated in AF-W (FNR=0.49) and AF-E (FNR=0.20).
116/117 GCH1 FNs have `segment_cn2_fraction=0` — the HMM never enters CN≥2 at those gene bins.
The FN CRR p50=1.42 shows the raw signal is present in coverage; the VAE reconstructs it, but not
anomalously enough to cross the HMM detection threshold.

Two targeted training changes should push the reconstruction error for weak amplifications above the
detection threshold:

1. **Boosted sinusoidal loss** (`sin_loss_max_weight` 1.0 → 2.0): the sinusoidal component
   penalises deviations from smooth periodic patterns more strongly. This makes the model more
   sensitive to block-level coverage changes like CNVs, increasing reconstruction error when a
   CNV is present but weak.

2. **Slower curriculum** (`cnv_downsample_warmup_epochs` 40 → 80): the weighted sampler ramps from
   5% to 30% CNV samples over 80 epochs instead of 40. The model spends more training time learning
   the normal coverage distribution before seeing CNV samples, tightening the normal-space boundary
   and making subtle amplifications like GCH1 in AF-W stand out more clearly as anomalies.

Both changes target the same root cause: the VAE's normal-space boundary is too permissive for weak
African amplifications. Switching to 1000bp bins (as originally proposed) would be counterproductive
— at 1000bp, GCH1 spans only ~2 bins, making the HMM call harder, not easier.

## Changes from experiment 36

| Parameter                      | Exp 36 | Exp 37   | Rationale                                                            |
|--------------------------------|--------|----------|----------------------------------------------------------------------|
| `sin_loss_max_weight`          | 1.0    | **2.0**  | Stronger sinusoidal penalty increases reconstruction error for weak CNVs |
| `cnv_downsample_warmup_epochs` | 40     | **80**   | Slower curriculum tightens normal-space boundary before CNV exposure |

All other parameters unchanged from exp 36. Full retrain (cannot reuse exp 32 checkpoint with
changed training regimen — training dynamics differ from epoch 1).

## Expected outcome

- **GCH1 FNR**: 0.16 → ≤0.10. Improved reconstruction sensitivity should push some AF-W/AF-E
  FNs above the HMM threshold. Not expecting full resolution in one step — the African amplification
  is genuinely weak.
- **CRT PPV**: ~0.92 (unchanged — calling parameters are identical)
- **MDR1 FNR**: should remain ~0.02
- **PM2_PM3 FNR**: should remain ~0.00
- **call_rate**: ~0.13 (noise filter unchanged; expect similar to exp 36)

## Run time estimate

~2–3 hours (full VAE training + HMM segmentation + CNV calling + evaluation).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 1.00 | 1.00 |
| GCH1    | 0.09 | 1.00 | 0.94 |
| MDR1    | 0.02 | 1.00 | 0.99 |
| PM2_PM3 | 0.00 | 1.00 | 1.00 |

**Predictions matched:** GCH1 FNR met the ≤0.10 target (0.09). CRT PPV stayed at 1.00 and MDR1 FNR
unchanged at 0.02, confirming that the training changes did not disturb other genes.

**AF-E resolved, AF-W improved:** AF-E GCH1 FNR dropped from 0.20 → 0.00 (fully resolved). AF-W
improved from 0.49 → 0.40, still the dominant failure mode. AS-S-FE also shows FNR=0.43 (n=16),
a smaller population with the same pattern.

**19 remaining GCH1 FNs — all frac=0, CRR p50=1.48:** The boosted sin loss raised FN CRR from
1.42 to 1.48 (the VAE now generates more reconstruction error for these weak amplifications).
But all 19 FNs still have segment_cn2_fraction=0 — the HMM never enters CN≥2, so the calling
threshold cannot recover them. CRR rescue (threshold 1.55) can't reach them either since p90=1.54.

**Key insight for exp 38:** The new VAE generates significantly lower CRR for TNs (GCH1 TN p90=1.10,
CRT TN p90=1.15, all genes ≤1.18). This creates ~17bp of headroom below 1.35 that didn't exist
with the old VAE. Lowering the rescue threshold to 1.35 is now safe — unlike exp 35, where old-VAE
CRT TNs reached CRR 1.52.

**Next experiment:** Exp 38 — lower `cnv_crr_rescue_threshold` 1.55 → 1.35 using exp 37 artifacts
(no retrain, no HMM re-run). Expected to rescue ~14–17 of the 19 remaining GCH1 FNs.

## Proposal history

**First proposal (2026-05-11):** Experiment 37 was originally proposed as a switch to 1000bp bins:
`store_path` 400bp → 1000bp, `hmm_self_transition` 0.90 → 0.80, `cnv_crr_rescue_threshold`
1.55 → 1.35, `cnv_max_transitions` 50 → 100. The rationale was that 1000bp bins aggregate 2.5x
more reads/bin and that historical 1000bp experiments achieved GCH1 FNR=0.04–0.07.

**Feedback received:** The 1000bp reasoning was flawed in two ways. First, the VAE is reconstructing
the GCH1 signal (CRR p50=1.42 is elevated) — the problem is that the reconstruction is not
anomalous enough, not that it fails. Second, at 1000bp, GCH1 spans only ~2 bins, making the HMM
call harder rather than easier. The right fix is to improve the VAE training regimen at 400bp —
specifically, boost the sinusoidal loss component and slow the curriculum warmup.

**Revised proposal:** Stay at 400bp. Two training changes: `sin_loss_max_weight` 1.0 → 2.0 and
`cnv_downsample_warmup_epochs` 40 → 80. Full retrain from scratch.

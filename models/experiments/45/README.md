# Experiment 45 — Flat 5% CNV curriculum: cnv_downsample_ratio_final 0.10 → 0.05

**Status:** Complete 2026-05-18

**Full retrain.** Curriculum change requires a new VAE from scratch; HMM segmentation and CNV
calling are also re-run.

## Hypothesis

Experiment 44 (10% CNV ratio) reduced GCH1 FNR from 0.15 → 0.09, recovering ~128 of the 306
FNs from exp 43. But 178 FNs remain, and 177/178 still have `segment_cn2_fraction=0`: the HMM
never assigns CN≥2 to any GCH1 bin for these samples. Their raw CRR p50=1.44 confirms real
signal exists, but the VAE reconstruction residual at GCH1 is still too low for the HMM to act.

The 10% ratio was not aggressive enough for a subpopulation of African GCH1-amplified profiles
(AF-E FNR=0.26, AF-W FNR=0.26): the VAE still saw enough amplified profiles to partially
learn them as near-normal. Setting the final ratio to 5% — equal to the initial ratio — eliminates
the curriculum ramp entirely. The VAE is exposed to just 5% CNV profiles throughout all 100 epochs,
maximising the reconstruction penalty for amplified profiles and giving the HMM a stronger signal.

## Changes from experiment 44

| Parameter                        | Exp 44 | Exp 45    | Rationale                                                                              |
|----------------------------------|--------|-----------|----------------------------------------------------------------------------------------|
| `cnv_downsample_ratio_final`     | 0.10   | **0.05**  | Flatten the curriculum at 5%; force VAE to treat amplification as genuinely anomalous  |
| VAE retrain                      | yes    | **yes**   | Curriculum change requires fresh training                                              |
| HMM re-run                       | yes    | **yes**   | New reconstructions require new segmentation                                           |

`cnv_downsample_ratio_initial` is already 0.05, so setting `final=0.05` eliminates the upward ramp.
All other parameters unchanged (β=4, `18_genome_cnv_caller`, global `cnv_crr_rescue_threshold=1.55`).

## Expected outcome

- **GCH1 FNR:** With a flat 5% curriculum, the VAE should reconstruct amplified profiles more
  poorly throughout training, generating larger reconstruction residuals at GCH1 bins and prompting
  the HMM to enter CN≥2 for some frac=0 FNs. If the mechanism holds, GCH1 FNR could drop
  toward 0.04–0.07.
- **HMM transitions p50:** May rise further from 209 (exp 44). The VAE will be more sensitive
  overall. Key diagnostic: if p50 exceeds ~350–400 without GCH1 FNR improvement, the 5% ratio
  is too aggressive and we should consider a different approach (excluded-normal training pool
  or noise filter re-enabled).
- **Callability:** Should remain near 0.996.
- **Other genes (CRT/MDR1/PM2_PM3):** FNR is already near-zero; small risk of PPV degradation
  if increased VAE sensitivity generates more false reconstruction residuals. Monitor MDR1 and
  PM2_PM3 PPV for regression.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.02 | 0.66 | 0.80 |
| GCH1    | 0.09 | 0.96 | 0.92 |
| MDR1    | 0.00 | 0.90 | 0.95 |
| PM2_PM3 | 0.02 | 0.86 | 0.91 |

HMM: callability=0.995, transitions p50=139 (down from 209 in exp 44).

**GCH1:** FNR unchanged at 0.09, but FN count increased marginally to 191 (vs 178 in exp 44). 190/191 FNs
have `segment_cn2_fraction=0` — the HMM never entered CN≥2 for any GCH1 bin in these samples. FN CRR
p50=1.44, confirming signal is present but the VAE reconstruction residual is insufficient.

**PPV dramatic improvement over exp 44:** CRT PPV 0.10→0.66, MDR1 0.51→0.90, PM2_PM3 0.53→0.86.
The flat 5% curriculum eliminated the ramp's excess sensitivity, reducing spurious HMM transitions.
This shows the flat curriculum was the right call for PPV — but it did not move GCH1 FNR.

**Divergence from expected:** GCH1 FNR was predicted to drop toward 0.04–0.07; it stayed at 0.09. The
transition p50 dropping (209→139) was predicted to rise — but this drop is actually a win (less
FP-generating over-segmentation). The core failure is that the VAE still reconstructs African
GCH1-amplified profiles with insufficient residual for the HMM to act, regardless of curriculum ratio.

**Conclusion:** Two consecutive curriculum adjustments (10% ramp in exp 44, flat 5% in exp 45) have
both plateaued at GCH1 FNR=0.09. With frac=0 for 190+ FNs, the evaluation framework confirms
parameter tuning alone won't fix this — the spatial-locality inductive bias of the conv encoder
must be removed. Next: experiment 46 switches to a fully-connected (MLP) VAE (07_fc_vae).

## Run time estimate

~7–8 hours: VAE training (~5–6h for 100 epochs) + HMM re-segmentation (~2h) + CNV calling
(~30 min) + evaluation (~5 min).

## Proposal history

**Original proposal (2026-05-17):** Reduce `cnv_downsample_ratio_final` from 0.10 → 0.05,
eliminating the curriculum ramp. Rationale: 177/178 remaining GCH1 FNs had `frac=0` in exp 44,
indicating the VAE still reconstructs a subset of African amplified profiles as near-normal.
Flat 5% was expected to maintain consistent reconstruction pressure throughout all 100 epochs.

**Feedback received:** Asked whether the 5% is sampled randomly each epoch, or always the same
fixed 5% of CNV samples.

**Answer:** Randomly sampled each epoch. `WeightedRandomSampler(weights, num_samples=len(dataset),
replacement=True)` assigns weight 0.05 to every CNV-positive sample and 1.0 to every normal
sample, then draws a fresh `len(dataset)`-sized index set with replacement at the start of each
epoch. Over 100 epochs, each CNV sample is seen ~5 times on average (vs ~100 times for a normal
sample) — a 20:1 frequency advantage for normal profiles — but the specific samples drawn vary
each epoch. No CNV sample is permanently excluded.

**What changed:** Nothing — the proposal stands. The answer confirms the mechanism: the VAE
eventually sees every CNV sample, just at 1/20th the frequency. The flat 5% is the correct
approach.

# Experiment 40 — Raise HMM self_transition 0.90 → 0.96 to fix sanity-check mass failure

**Status:** Complete 2026-05-13

**No retrain, no inference re-run.** Reuses exp 37 checkpoint, latents, and reconstructions
via symlinks. HMM segmentation is re-run with the new self_transition. CNV calling and
evaluation are re-run.

## Hypothesis

Experiment 39 removed the noise filter and revealed a new blocker: **the sanity check**.
47–67% of ground-truth CNV positives are still called -1. They have high CRR (GCH1 masked
p50=2.11, all genes p25 > 1.55), confirming the amplification signal is present. The
sanity check requires ≥55% of chromosome bins to be covered by high-confidence CN=1
segments; with transitions p50=1189, the HMM oscillates so rapidly that many chromosomes
fall below 55% CN=1, causing the sanity gate to fire even on true amplification samples.

**Root cause:** `hmm_self_transition=0.90` is too low for the current 06_conv_vae model.
At exp 30 with an older VAE the same 0.90 produced p50=270 transitions; the new model
produces noisier reconstruction residuals that drive the HMM to switch states much more
frequently, giving p50=1189. The sanity check is correct in principle — it is the HMM
calibration that is wrong.

**Proposed fix:** Raise `hmm_self_transition` from 0.90 to 0.96. This value was used
in exp 30/31 and produced clean segmentation. The reason it was later abandoned:
exp 31 showed high evaluable FNR (GCH1=0.49) because with a sticky HMM, amplifications
were not segmented (segment_cn2_fraction=0 for most FNs) and there was **no CRR rescue**
at that time. CRR rescue was added in exp 33/35. Now, even if the stickier HMM fails to
segment an amplification body, the rescue criterion (CRR ≥ 1.55) catches it directly
from raw coverage.

## Changes from experiment 39

| Parameter              | Exp 39 | Exp 40   | Rationale                                                                    |
|------------------------|--------|----------|------------------------------------------------------------------------------|
| `hmm_self_transition`  | 0.90   | **0.96** | Reduces over-segmentation; fixes sanity-check mass failure                   |
| HMM re-run             | no     | **yes**  | Self-transition change requires new segments.parquet                         |

All other parameters unchanged. Noise filter remains disabled (`cnv_max_transitions: null`).

## Expected outcome

- **HMM transitions p50**: drops from 1189 toward O(10–100), as in exp 30.
- **Sanity check pass rate**: most samples pass (chromosome CN=1 proportion returns to
  normal), dramatically reducing the -1 pool.
- **masked_cnv_rate**: drops from 47–67% toward a small residual. Samples that genuinely
  have low coverage or very noisy chromosomes will still fail sanity, but should be rare.
- **CRR rescue role**: for CNV positives where the stickier HMM produces segment_cn2_fraction=0,
  CRR rescue fires instead — masked CRR p25=1.61–2.00 across all genes ensures most are
  above the 1.55 threshold.
- **total_fnr**: should drop from 0.40–0.71 toward a value much closer to the evaluable
  FNR (0.02–0.11).
- **PPV**: the evaluable pool grows as more samples pass sanity; PPV may shift.
  FP risk is limited — for sanity-passing samples, TN CRR p90=1.07–1.24 (well below 1.55).
- **CRT**: likely still challenging — few GT positives (58 total) and high FP rate.

## Run time estimate

~2–3 hours: HMM re-segmentation (~2h with n_jobs=-1 across 50k samples) +
CNV calling (~30 min) + evaluation (~5 min). No training or inference.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.11 | 0.02 | 0.11 |
| GCH1    | 0.09 | 0.59 | 0.69 |
| MDR1    | 0.02 | 0.57 | 0.74 |
| PM2_PM3 | 0.04 | 0.37 | 0.57 |

Masked CNV rate: CRT=0.67, GCH1=0.47, MDR1=0.40, PM2_PM3=0.52.
Total FNR: CRT=0.71, GCH1=0.52, MDR1=0.41, PM2_PM3=0.54.
HMM transitions p50=970 (predicted O(10–100); actual only 18% reduction from 1189 in exp 39).

**Where predictions matched:** Evaluable FNR is excellent (0.02–0.11), confirming CRR rescue
is effective for the subset that passes the sanity check. The sanity check blocker diagnosis
was correct.

**Where predictions diverged:** The key prediction was that 0.96 would restore clean
segmentation as it did in exp 30/31 with an older VAE. It did not. The 06_conv_vae residuals
are fundamentally noisier than the old model: at the same self_transition of 0.96, the new VAE
yields p50=970 vs the old VAE's p50=270. The transition probability is 4× lower (0.04 vs 0.10)
but transitions only dropped ~18% — the HMM emission variance from the new VAE is large enough
that the state-transition prior barely constrains oscillation. The sanity check (≥55% CN=1)
continues to fire on the majority of both CNV-positive and normal samples with complex coverage.

**Next:** Experiment 41 — raise self_transition further to 0.99 to target p50~240 transitions.

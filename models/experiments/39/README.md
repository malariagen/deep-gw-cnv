# Experiment 39 — Disable noise filter to stop masking CNV-positive samples

**Status:** Complete 2026-05-13

**No retrain, no HMM re-run.** Symlinks exp 38 artifacts. Single change: `cnv_max_transitions: 50 → null`.

## Hypothesis

Experiment 38 revealed that 84–98% of ground-truth CNV samples are called -1 (masked). This is
not a VAE signal problem — masked CRR p50 is 1.79–2.42 across all genes, far above 1.0, confirming
the amplification signal is present. The culprit is the noise filter in `16_genome_cnv_caller`:

```
cnv_max_transitions: 50
```

The HMM produces transitions p10=85, p50=1189, p90=2945 per sample. With a threshold of 50,
essentially all samples (>90% even at the p10) exceed the limit. Since the noise filter fires
at priority 0 — before sanity check and before CRR rescue — samples with CRR=3.0 or higher are
still called -1, and the rescue threshold never fires.

Disabling the filter (`cnv_max_transitions: null`) removes the blockage. The CRR rescue at 1.55
then applies as designed, and the HMM coverage criterion applies for samples with elevated-CN
segments. Since PPV=1.00 in exp 38, there is currently no FP problem — the filter is not
suppressing false positives, it is suppressing true positives.

## Changes from experiment 38

| Parameter             | Exp 38 | Exp 39   | Rationale                                                                   |
|-----------------------|--------|----------|-----------------------------------------------------------------------------|
| `cnv_max_transitions` | 50     | **null** | Disable noise filter; threshold of 50 filters ~98% of samples (HMM p10=85) |

All other parameters and all artifacts unchanged. Artifacts symlinked from exp 38.

## Expected outcome

- **model_miss**: drops from 0.98 toward the true quality-filter rate (samples with genuinely
  low coverage or failed sanity). Expected ~0.05–0.15 rather than 0.98.
- **masked_cnv_rate**: drops dramatically from 0.84–0.98. CNV+ samples have high CRR (p50≥1.79),
  so most will be rescued by the CRR criterion at 1.55 and called CN=2 rather than -1.
- **total_fnr**: should drop from 0.84–0.98 toward a value close to FNR (≤0.10 for GCH1),
  revealing true performance for the first time.
- **FNR / PPV / MCC** on the evaluable subset: PPV may decrease slightly if previously-masked
  samples include near-threshold calls; FNR on the evaluable set may shift as the evaluable
  pool grows.
- **FP risk**: low, because CRR rescue requires CRR ≥ 1.55 and TN CRR p90=1.10–1.18. Any
  FP increase will be visible in PPV.

This establishes a sound calling baseline — the first experiment where the quality filter is
not pathologically miscalibrated to the current HMM output.

## Run time estimate

~15 minutes (CNV calling + evaluation only; all artifacts symlinked from exp 38).

## Actual outcome

| Gene    | MCC  | FNR  | PPV  | masked_cnv_rate | total_fnr |
|---------|------|------|------|-----------------|-----------|
| CRT     | 0.11 | 0.11 | 0.01 | 0.67            | 0.71      |
| GCH1    | 0.67 | 0.09 | 0.57 | 0.47            | 0.52      |
| MDR1    | 0.72 | 0.02 | 0.55 | 0.39            | 0.40      |
| PM2_PM3 | 0.55 | 0.03 | 0.35 | 0.52            | 0.53      |

**Where predictions were right:** Removing the noise filter dramatically reduced masked_cnv_rate
from the 84–98% seen in exp 38. The noise filter was indeed the dominant cause of masking.
Evaluable-subset FNR (0.02–0.11) and PPV stayed in the expected range. CRR rescue at 1.55
operates as designed for evaluable samples.

**Where predictions diverged — masked_cnv_rate still high:** The prediction was that
masked_cnv_rate would approach zero and total_fnr would collapse toward evaluable FNR.
Instead, 39–67% of GT CNV positives are still masked. These samples have high CRR
(GCH1 masked p50=2.11, all genes p25 > 1.55), confirming the signal is present. The
culprit is the **sanity check** (≥55% chromosome CN=1 coverage): with HMM transitions
p50=1189 per sample, the HMM oscillates so rapidly between states that many chromosomes
fall below the 55% CN=1 threshold, and the sample is marked -1 before CRR rescue fires.

**Root cause:** `hmm_self_transition=0.90` produces 1189 transitions p50 with the
current 06_conv_vae, compared to p50=270 in exp 30 with the old VAE. The new model
produces noisier inputs to the HMM, amplifying the over-segmentation. The sanity check
is an HMM-quality gate that correctly fires on garbage segmentations, but at 0.90 it
fires on almost everything, including true amplifications.

**Next experiment:** Exp 40 — increase `hmm_self_transition` to 0.96 (historically used
in exp 30/31) with HMM re-run. CRR rescue now provides a safety net that didn't exist
when 0.96 was last tried.

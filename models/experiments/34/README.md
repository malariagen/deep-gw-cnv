# Experiment 34 — Per-sample noise filter + CRR rescue (v16) on exp 32 artifacts

**Status:** Complete 2026-05-06

**No retraining, no new HMM run.** Symlinks exp 32 checkpoint, reconstructions, and
segments. Changes are entirely in the CNV caller (v16).

## Hypothesis

Experiments 31–33 converge on the same finding: 70–88% of FNs have
`segment_cn2_fraction=0` regardless of HMM parameters or training depth. The CRR signal
is clearly present (FN CRR p50 = 2.24–2.66) but the HMM cannot segment gene-scale
amplifications.

Before applying CRR rescue, however, we must address a more fundamental problem:
the HMM is over-segmenting aggressively. Exp 32 shows transitions p50=270 (per sample,
summed across all chromosomes). A haploid organism with a few CNVs at a time should
have O(10) CN state changes; hundreds of transitions indicate the HMM is fitting noise,
not signal. Calling on these samples inflates FPs and makes the CRR rescue threshold
less discriminative (noisy samples may also have locally elevated CRR from HMM
fragmentation artifacts).

The v16 caller adds two changes:
1. **Noise filter**: before any criterion fires, samples with total CN transitions > 50
   (roughly the 15th percentile of exp 32) are marked CN=-1. This restricts calling to
   samples where the HMM produced a coherent, low-fragmentation segmentation.
2. **CRR rescue** (retained from original v15 proposal): for clean samples that still
   fail the HMM criterion, rescue if gene/flank CRR ≥ 1.60. FP CRR p90 = 1.50–1.53,
   FN CRR p25 = 1.67–2.23 — the threshold sits cleanly between them.

## Changes from experiment 33

| Component              | Exp 33                          | Exp 34                                    | Rationale                                                    |
|------------------------|---------------------------------|-------------------------------------------|--------------------------------------------------------------|
| `cnv`                  | `14_genome_cnv_caller`          | **`16_genome_cnv_caller`**                | Adds noise filter + CRR rescue                               |
| `cnv_max_transitions`  | N/A                             | **50**                                    | Filter samples where HMM is fitting noise (~85% excluded)    |
| `cnv_crr_rescue_threshold` | N/A                         | **1.60**                                  | Above FP p90 (1.53); below FN p25 (1.67–2.23) for clean samples |
| Checkpoint/recons      | Retrained (200 ep)              | **Symlinked from exp 32**                 | Exp 33 retraining made results worse; revert                 |
| `segments.parquet`     | Exp 33 HMM (transitions p50=406)| **Symlinked from exp 32**                 | Exp 32 HMM less fragmented (p50=270), better PPV baseline    |
| `out_dir`              | `..._200ep_st0.90`              | **`..._st0.90_noise_filter_crr_rescue`**  | Reflects both changes                                        |

## Expected outcome

- **Callable sample rate**: ~15% of samples pass the noise filter (transitions ≤ 50).
  Missingness will increase substantially. This is by design: we are trading coverage
  for call quality on the samples that are called.
- **PPV**: should improve or hold. Filtering out noisy samples removes the main source
  of HMM-fragmentation FPs. CRR rescue on clean samples has a cleaner baseline
  (lower TN CRR p90 expected after noise filtering), so the 1.60 threshold is safer.
- **FNR**: unclear. FNs in the clean subset may still have segment_cn2_fraction=0 but
  high CRR — the rescue criterion should recover those. The 85% filtered out will
  all be CN=-1 (uncallable), which the evaluation treats as missed. Whether this is
  better than the current outcome depends on how many TPs are in the noisy vs clean subset.
- **MCC**: net improvement expected if the filtered sample set has higher TNR and the
  rescue recovers enough TPs.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.04 | 0.92 | 0.94 |
| GCH1    | 0.20 | 0.99 | 0.88 |
| MDR1    | 0.02 | 0.89 | 0.93 |
| PM2_PM3 | 0.00 | 0.90 | 0.95 |

**Missingness delta**: +0.29–0.32 (model_miss ≈ 0.87 vs pf9_miss ≈ 0.57). Consistent with
the ~85% noise filter prediction; actual call_rate ≈ 13%.

**What matched predictions:** PPV improved dramatically (0.89–0.99) vs exp 32 (0.53–0.70).
CRR rescue is working — 593 GCH1 TPs, 317 MDR1 TPs, 254 PM2_PM3 TPs caught by rescue.
TN CRR background for callable samples is tight (GCH1 p90=1.12).

**Where results diverged:** GCH1 FNR=0.20 (147 FNs) is higher than expected. Near-miss
diagnostics show 146/147 FNs have segment_cn2_fraction=0 — HMM gave no coverage. But
CRR BY CALL OUTCOME shows FN p50=1.48, p90=1.58: these FNs have genuine CRR signal, just
below the 1.60 rescue threshold. The threshold was calibrated on exp 32's noisy TN
background (TN p90=1.50–1.53); in the cleaner callable subset, TN p90=1.12, so 1.60 is
now overly conservative. Lowering to 1.45 would recover most GCH1 FNs.

**Next:** Experiment 35 — lower `cnv_crr_rescue_threshold` from 1.60 to 1.45.

## Run time estimate

~30 minutes (CNV calling + evaluation only; all artifacts symlinked from exp 32).

## Proposal history

**Original proposal (2026-05-06):** CRR rescue caller (v15) on exp 32 artifacts, with
`cnv_crr_rescue_threshold: 1.60`. Rationale: exps 31–33 show HMM cannot segment gene-scale
amplifications; FN CRR p50 = 2.24–2.66 clearly separable from FP CRR p90 = 1.50–1.53.

**Feedback 1 (2026-05-06):** Terminology correction — "diploid" → "haploid". Experiment
design unchanged.

**Feedback 2 (2026-05-06):** A sample shouldn't have more than a few CNVs at a time. If
there are too many transitions, the sample is too noisy to call. Stop trying to perceive
everything as signal. TN CRR should be as close to 1 as possible; TP CRR should be close
to 2 or higher.

**Revision (2026-05-06):** Added per-sample noise filter (`cnv_max_transitions: 50`).
New cnv version v16 (cannot edit v15 in place per audit-trail rules). Rationale: exp 32
transitions p50=270, p90=1796 — the vast majority of samples are too fragmented. Filtering
at 50 (≈p15) keeps only samples with coherent HMM segmentation before applying CRR rescue,
so the rescue fires only where TN CRR background is actually low. `out_dir` updated to
reflect both changes.

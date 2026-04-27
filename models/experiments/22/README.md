# Experiment 22 — Synthetic normal augmentation (Poisson resampling)

**Status:** Complete 2026-04-21

**Full retraining required: VAE fine-tuned from exp 21 checkpoint with Poisson-resampled normal samples.**

## Hypothesis

Experiments 19–21 leave GCH1 FNR stuck at 0.07 with FN CRR p50=1.25. The signal is genuine
(p50 >> 1.0), but the VAE reconstructs CNV profiles partially faithfully because it has seen
CNV samples during training. Lowering the gate threshold to catch CRR=1.25 is not the right
fix — it chases weak signal and trades FPs for FNs (CRR=1.25 is uncomfortably close to the
normal baseline of 1.0).

The correct fix is to boost the signal: force the VAE to produce higher reconstruction error
for true CNV samples so their CRR rises above 1.30 on its own.

**Mechanism**: During training, every normal (non-CNV) sample is Poisson-resampled in
`__getitem__` (fresh draw each epoch). This gives the VAE a much richer picture of what
normal coverage variation looks like. As a result:
1. The VAE encodes normal samples into a tighter, more precise region of latent space.
2. At inference, true CNV samples fall further from any normal archetype → decoder
   outputs a more "normal-looking" reconstruction → reconstruction error (counts/recon)
   is larger → CRR increases.

CNV-positive samples are passed through unchanged; curriculum sampling (ratio 0.05→0.30)
is retained so the model still sees enough CNV examples to segment correctly.

## Changes from experiment 21

| Component                       | Exp 21                    | Exp 22                            | Rationale                                                  |
|---------------------------------|---------------------------|-----------------------------------|------------------------------------------------------------|
| `aug_normal_poisson`            | not set                   | `true`                            | Poisson-resample normals each epoch; tighter normal prior  |
| `aug_normal_depth_scale`        | not set                   | `[0.85, 1.15]`                    | Also randomly scale depth; adds realistic library-size variation |
| VAE checkpoint (start)          | exp 18 (fine-tune from)   | exp 21 (fine-tune from)           | Warm start from most recent good checkpoint                |
| `cnv_crr_amp_threshold`         | 1.30                      | 1.30 (unchanged)                  | Gate stays put; signal-boosting via training, not gate     |
| Everything else                 | unchanged                 | —                                 | Same architecture, sin loss, curriculum, HMM params        |

## Expected outcome

- **GCH1 FNR**: if augmentation succeeds in tightening the normal prior, FN CRR should shift
  upward. Target: FNR from 0.07 toward 0.03–0.05. AF-E (FNR=0.27) and AF-W (FNR=0.18) should
  benefit most if their CRR is genuinely in [1.25, 1.30).
- **GCH1 PPV**: augmentation adds normal training variation; the model should not inflate
  normal-sample CRR (those samples are seen in their varied form, not as CNVs). Risk: if
  Poisson noise pushes some normal samples' CRR above 1.30, PPV may decrease slightly.
  TN p90=1.13 provides a buffer; moderate Poisson noise should stay well clear of 1.30.
- **PM2_PM3 / CRT / MDR1**: no expected material change; curriculum ratio unchanged.

## Run time estimate

~2–3 hours (VAE fine-tuning from exp 21 checkpoint + HMM + CNV calling + evaluation on Mac mini).

## What we could do instead

1. **Masked reconstruction loss on gene regions** — zero the VAE reconstruction loss on bins
   overlapping known CNV gene regions during training, so the model is not penalised for failing
   to reconstruct amplifications. Requires a code change to compute_loss; stronger mechanism.
2. **Higher depth-scale range (0.70–1.30)** — more aggressive depth variation could further
   tighten the normal prior, but risks destabilising training.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| GCH1    | 0.06 | 0.87 | 0.90 |
| CRT     | 0.00 | 0.97 | 0.98 |
| MDR1    | 0.01 | 0.97 | 0.94 |
| PM2_PM3 | 0.00 | 0.94 | 0.90 |

GCH1 FN CRR p50=1.23 (116 FNs). AF-E FNR=0.32, AF-W FNR=0.14.

**Where prediction matched:** CRT/MDR1/PM2_PM3 unchanged as expected. PPV risk was correctly
flagged — augmentation did add some noise to normal samples, inflating CRR slightly (PPV dropped
0.91→0.87).

**Where prediction diverged:** Augmentation did not boost FN CRR. FN p50 actually decreased
(1.25→1.23) rather than rising above 1.30. AF-E FNR got *worse* (0.27→0.32), the opposite of
the predicted improvement. The Poisson resampling tightened the normal side of the distribution
slightly (fewer TNs above 1.15) but had no meaningful effect on the CNV side.

**Why it diverged:** The CNV signal captured in the latent space appears robust to normal-sample
augmentation diversity. The model already had a sufficiently tight normal prior — more diverse
normals during training didn't further separate CNV from normal in reconstruction space. The FNs
at CRR=1.23 are genuinely hard cases that look similar to normals; augmentation cannot fix that.

**Next:** Experiment 23 — genome-wide CNV calling across all 5,318 P. falciparum protein-coding
genes using `06_genome_cnv_caller`. Threshold lowered to 1.25 (same rationale as original exp 23
proposal). Produces a master amplification table (one row per amplified sample × gene) in addition
to the standard four-gene evaluation.

## Proposal history

**Original proposal (2026-04-20, before feedback):** Lower `cnv_crr_amp_threshold` from 1.30
to 1.25 to capture FNs with CRR in [1.25, 1.30). No retraining; symlinks from exp 21.

**Feedback received:** CRR=1.25 is too close to the normal baseline (1.0) to safely lower the
gate. The underlying problem is that the VAE signal is weak, not that the gate is too high.
Training-level signal boosting (synthetic augmentation of normal samples) was specifically
requested as the correct approach.

**What changed:** Completely redesigned from a gate-lowering experiment to a training experiment.
`cnv_crr_amp_threshold` stays at 1.30. New `aug_normal_poisson=true` Poisson-resamples each
normal sample per-epoch during fine-tuning. Full retraining from exp 21 checkpoint required.

**Crash and fix (2026-04-15):** Experiment was authorised and started but crashed immediately
with `AttributeError: 'AugmentedNormalDataset' object has no attribute 'sample_ids'`. The
curriculum loader in `train.py` accesses `ds.sample_ids` after `ds` has been wrapped in
`AugmentedNormalDataset`. Fixed by adding `self.sample_ids = base_dataset.sample_ids` to
`AugmentedNormalDataset.__init__` in `models/training/dataset.py`. No config or design changes.

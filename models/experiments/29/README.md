# Experiment 29 — 06_conv_vae on 400 bp core genome (dynamic n_bins)

**Status:** Proposed 2026-04-27

**Full training from scratch.** Architecture 06_conv_vae with 400 bp bin resolution
(51986 bins). No checkpoint reuse — the linear layer dimensions change with bin count.

## Hypothesis

Experiments 21–28 all train or infer on 1000 bp bins (20814 bins covering the core
genome). At 1000 bp, a single bin spans a kilobase — short-gene amplifications (GCH1
is only ~1.9 kb, about 2 bins) produce a very coarse anomaly signal. The HMM then sees
an elevation over just 1–2 positions, which makes partial-coverage FNs more likely.

Switching to 400 bp bins gives ~2.5× more positional resolution (51986 bins). GCH1 now
spans ~5 bins, CRT ~10, MDR1 ~18. The reconstruction error is more spatially precise,
and the HMM has more consecutive elevated-CRR bins to commit to. This should improve
callability for short amplifications without touching the HMM or caller parameters.

The architecture change (06_conv_vae) is purely mechanical: `ConvVAE(n_bins, latent_dim)`
receives `n_bins` from the dataset at runtime instead of reading a module-level constant.
The topology (5× ResConvBlock encoder + 5× ConvTranspose1d decoder) is unchanged from 05.

## Changes from experiment 28 (proposed) / 27 (last completed)

| Component             | Exp 27/28 | Exp 29 | Rationale                                               |
|-----------------------|-----------|--------|---------------------------------------------------------|
| `architecture`        | 05_conv_vae | 06_conv_vae | n_bins auto-derived from data; same topology      |
| `store_path`          | 1000 bp core | 400 bp core | 2.5× more positional resolution                  |
| `batch_size`          | 128       | 64     | Halved: wider inputs increase per-batch memory          |
| Training              | no (reuse checkpoint) | full | Different bin count → linear layers differ    |
| HMM / CNV caller      | same      | same   | CRR is resolution-agnostic; carry over thresholds       |

## Infrastructure changes (one-off, shared across future experiments)

- `architectures/06_conv_vae.py`: n_bins is a constructor param, not a constant.
- `training/dataset.py`: `ReadCountDataset` derives `n_bins_raw` / `n_bins_padded`
  from data shape; loads counts with `mmap_mode='r'` (saves ~11 GB RAM at init).
- `train.py` / `wrap_up.py`: pass `n_bins=ds.n_bins_raw` to architectures that support it.

## Expected outcome

- **GCH1 FNR**: 0.11 → likely 0.04–0.08. More bins per gene means HMM sees a longer
  elevated run; marginal-CRR TPs (1.25–1.35) should have stronger CN=2 segment coverage.
- **MDR1 FNR**: 0.01 → ~0.01. Already low; no regression expected.
- **PM2_PM3 FNR**: 0.06 → unclear. The PM2/PM3 FNs at CRR≈0 are coverage failures,
  not resolution failures; no improvement expected there. Short-span TPs may improve.
- **PPV**: Starting fresh means the new model may have different systematic biases.
  The two-tier band caller is still in place so artefact suppression should hold;
  expect PPV in a similar range to exp 27.

## Run time estimate

~6–10 hours (full training: 100 epochs × 53973 samples × 51986 bins at batch=64,
plus HMM + CNV calling + evaluation). HMM with n_jobs=-1 will be ~2.5× slower than
1000 bp due to more positions per chromosome.

## Proposal history

**Original proposal (2026-04-27):** Proposed switching from 1000 bp to 400 bp bins with
architecture 06_conv_vae (n_bins auto-derived), halved batch size, and full retraining
from scratch. Two-tier band caller (v11) carried over from exp 28 proposal. Motivated
by GCH1 FNR regression to 0.11 (from 0.04) in exp 27: at 1000 bp, GCH1 spans only
~2 bins, giving the HMM an imprecise anomaly signal.

**Feedback received (2026-04-27):** Two questions — (1) explanation of the band filter
requested; (2) confirm that data augmentation is not being used.

**Revised proposal (2026-04-27):** No changes to config or architecture. Added clearer
explanation of the two-tier band filter concept in the accompanying email. Confirmed
that `AugmentedNormalDataset` (Poisson resampling) is not active — exp 29's config
does not set `aug_normal_poisson`, so training runs on raw log2-normalised read counts.

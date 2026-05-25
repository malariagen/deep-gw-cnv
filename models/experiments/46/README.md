# Experiment 46 — Fully-connected (MLP) VAE: remove spatial-locality bias

**Status:** Re-proposed 2026-05-21 (bugfix: shape mismatch in 07_fc_vae.py)

**Full retrain.** Architecture change requires fresh training; VAE, HMM, and CNV calling all re-run.

## Hypothesis

Experiments 44 and 45 both plateau at GCH1 FNR=0.09 with 190/191 FNs having
`segment_cn2_fraction=0`. The VAE reconstruction residual is simply not large enough for
these samples — no amount of HMM or threshold tuning will fix it.

The root cause is that prevalent African GCH1 amplifications (AF-E FNR=0.36, AF-W FNR=0.26)
contaminate the 95% "normal" training pool. The convolutional VAE can learn "elevated counts
around a specific genomic region is normal in African samples" because its local receptive
fields can specialise on GCH1 bins without disturbing the rest of the genome. The network
effectively memorises the spatial pattern.

A fully-connected VAE applies the same global weight matrix to the entire coverage vector.
There is no local receptive field — the network cannot attend to GCH1 bins in isolation.
Any reconstruction must be mediated through the same latent space that encodes the global
coverage shape. This should make it structurally harder to represent population-specific
regional amplification patterns as "normal", raising reconstruction residuals for GCH1-amplified
profiles even at β=4.

## Changes from experiment 45

| Parameter     | Exp 45          | Exp 46              | Rationale                                              |
|---------------|-----------------|---------------------|--------------------------------------------------------|
| `architecture`| `06_conv_vae`   | **`07_fc_vae`**     | Remove spatial-locality inductive bias                 |
| `max_beta`    | 4.0             | 4.0 (unchanged)     | Isolate architecture as the single variable            |
| `out_dir`     | `45_06convvae…` | **`46_07fcvae…`**   | Reflect revised architecture                           |

Flat 5% curriculum retained from exp 45: it gave dramatic PPV improvements (CRT 0.10→0.66,
MDR1 0.51→0.90, PM2_PM3 0.53→0.86) and keeping it here means architecture is the only
change.

## Expected outcome

- **GCH1 FNR:** With no local receptive field, the VAE should struggle to reconstruct
  African GCH1-amplified profiles through the global latent space, raising reconstruction
  residuals and prompting the HMM to enter CN≥2 for some of the 190 frac=0 FNs.
  If the mechanism holds, GCH1 FNR could drop toward 0.05–0.07.
- **PPV for other genes:** Should remain near exp 45 levels (CRT=0.66, MDR1=0.90,
  PM2_PM3=0.86). Risk: if the FC encoder fits the normal coverage poorly overall
  (too noisy a reconstruction), transitions p50 will rise and PPV will drop.
- **HMM transitions p50:** Was 139 in exp 45. Monitor: if it rises above ~400 without
  GCH1 FNR improvement, the FC architecture is producing noisy reconstructions for
  normal samples (not useful signal).
- **Training time:** Likely faster than the conv VAE: the FC encoder/decoder is
  computationally simpler per forward pass, though the first layer (51986 × 2048)
  is large.

## What we could do instead

1. **FC VAE + β=8:** If the FC architecture alone is insufficient, stack with the
   originally proposed β=8 to add both an architectural and KL regularisation change.
2. **Exclude GT-positive samples from normal pool:** Use Pf9 ground truth to remove
   known GCH1-amplified samples from the 95% "normal" curriculum pool. Directly
   addresses the contamination problem but requires training data loader changes.

## Run time estimate

~6–8 hours: VAE training (~5–6h for 100 epochs) + HMM re-segmentation (~2h) + CNV calling
(~30 min) + evaluation (~5 min).

## Proposal history

**Original proposal (2026-05-18):**
Raise `max_beta` from 4.0 to 8.0. The argument was that stronger KL regularisation
forces the latent space toward a tighter global normal distribution, preventing the
VAE from dedicating latent capacity to population-specific patterns like prevalent
African GCH1 amplifications. Architecture unchanged from exp 45.

**Feedback received (2026-05-19):**
Chiyun asked whether the project could explore fully connected layers instead.

**What changed and why:**
The architectural change is a more fundamental intervention than a β increase: it
removes the spatial locality inductive bias entirely rather than just tightening the
KL regularisation. A FC VAE cannot rely on local filters to memorise regional patterns,
addressing the root cause more directly. β was reset to 4.0 to isolate the architecture
as the sole variable.

**Bug and fix (2026-05-21):**
The first run attempt failed immediately with `RuntimeError: linear(): input and weight.T
shapes cannot be multiplied (64x52000 and 51986x2048)`. `ReadCountDataset` pads every
sample to `ceil(n_bins/32)*32 = 52000` before the DataLoader, but `07_fc_vae.py`'s first
`Linear` layer was built with `n_bins_raw=51986`. Fixed in `07_fc_vae.py` by computing
`n_bins_padded` internally (same pattern as `06_conv_vae.py`) and using it for all linear
layer dimensions. `train.py` is unchanged — it still passes `n_bins_raw` and each
architecture handles padding internally.

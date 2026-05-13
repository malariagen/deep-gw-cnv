# Experiment 33 — Deeper VAE training (200 epochs) to improve reconstruction quality

**Status:** Complete 2026-05-05

**Full retraining from scratch.** Increases training epochs from 100 to 200 (patience 20 → 30).
Full pipeline: training → inference → HMM (st=0.90) → CNV calling (v14) → evaluation.

## Hypothesis

Experiments 29–32 all reuse the same checkpoint trained to 100 epochs. The VAE reconstruction
is systematically 0.5–0.8× the input signal — the model hit the epoch cap before fully
converging. This underprediction inflates CRR = coverage / reconstruction uniformly across all
samples, compressing the contrast between amplified and normal regions that the HMM relies on.

If the VAE converges further (reconstruction closer to input for diploid regions), then:
- CRR ≈ 1.0 for normal genomic regions in non-amplified samples
- CRR >> 1.0 specifically at amplified gene loci in CNV-positive samples

This cleaner signal should make it easier for the HMM to detect the gene-scale copy-ratio
elevations that exp 32 confirmed are present in the data (FN CRR p50: GCH1=1.84, MDR1=2.66,
PM2_PM3=2.53) but are not being segmented.

## Changes from experiment 32

| Component              | Exp 32                         | Exp 33                         | Rationale                                              |
|------------------------|--------------------------------|--------------------------------|--------------------------------------------------------|
| `epochs`               | 100 (hard cap, reused ckpt)    | **200**                        | Allow full convergence; reconstruction was 0.5-0.8× input |
| `patience`             | 20                             | **30**                         | Avoid premature stopping during slow late-phase improvement |
| `cnv`                  | `14_genome_cnv_caller`         | same                           | No rescue needed — testing reconstruction improvement alone |
| Checkpoint             | Reused from exp 29             | **Retrained from scratch**     | New training run required                              |
| HMM                    | Re-run at st=0.90              | **Re-run at st=0.90**          | Same setting; full re-run required after new checkpoint |
| `out_dir`              | `..._st0.90_hmm_primary_conf0.25` | **`..._200ep_st0.90`**     | Reflects the deeper training                           |

## Expected outcome

- **Reconstruction quality**: output/input ratio improves toward 1.0 for diploid regions.
  CRR p90 for normal (TN) samples should drop from 1.19–1.24 toward 1.05–1.10, sharpening
  the contrast with CNV-positive FN CRR values (1.84–2.66).
- **FNR**: if the HMM can now segment gene-scale events it previously missed, FNR should
  improve substantially for MDR1 and PM2_PM3 (where FN CRR >> 1.60). GCH1 remains uncertain
  (only 4–5 bins at 400bp; gene-scale fitting problem may persist).
- **PPV**: should be stable or improve slightly — better reconstruction reduces uniformly
  elevated CRR that could cause false CN>=2 calls.
- **MCC**: net improvement expected if FNR recovery is meaningful.

## Run time estimate

~8–10 hours (training ~7–8h; HMM ~2h; calling + evaluation ~30 min; steps overlap in
wrap_up.py pipeline triggered by train.py).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.25 | 0.05 | 0.19 |
| GCH1    | 0.47 | 0.48 | 0.44 |
| MDR1    | 0.27 | 0.23 | 0.37 |
| PM2_PM3 | 0.40 | 0.23 | 0.33 |

Deeper training made things worse across the board compared to exp 32 (FNR flat or
slightly higher; PPV collapsed from 0.14–0.70 to 0.05–0.48; MCC down significantly).

The reconstruction hypothesis was not confirmed: TN CRR p90 is still 1.18–1.25,
essentially identical to exp 32. The VAE did not produce a more accurate diploid
reconstruction despite 200 epochs. What changed is the HMM behaviour: transitions p50
increased from 270 (exp 32) to 406 (exp 33), meaning the new checkpoint causes the HMM
to over-segment more aggressively on genome-wide residuals, generating more FPs without
fixing detection at gene loci.

Near-miss frac=0 rates remain 70–88% across genes — the fundamental problem is unchanged.
FN CRR p50 values (CRT=2.44, GCH1=2.24, MDR1=2.59, PM2_PM3=2.66) are still well above
FP CRR p90 (1.50–1.53), confirming the signal is there but the HMM cannot segment it
regardless of training depth.

→ Experiment 34 (revised): CRR rescue (v16) on exp 32 artifacts, with additional per-sample
  noise filter (cnv_max_transitions=50). Feedback on original exp 34 proposal highlighted that
  samples with many HMM transitions are too noisy to call and should be excluded rather than
  rescued.

## Proposal history

**Original proposal (2026-04-29):** Exp 33 was proposed as a CRR rescue fallback — no
retraining, using v15 caller with `cnv_crr_rescue_threshold=1.60` to call CN=2 when HMM
misses but gene/flank CRR ≥ 1.60. Estimated 30 minutes (caller + evaluation only; segments
symlinked from exp 32).

**Feedback received:** "Can we train deeper? We finish at 100 epochs, and I feel the
reconstruction baseline signal is usually a factor of 0.5-0.8 off of the input."

**Revised approach:** The feedback identifies a more upstream problem: the VAE reconstruction
itself is poor, not just the HMM downstream. If the reconstruction underpredicts by 20-50%
uniformly, the CRR contrast between normal and amplified regions is compressed, making the
HMM's segmentation task harder. Fixing reconstruction quality first is the cleaner
intervention — it avoids adding a bypass (CRR rescue) on top of a flawed representation.
The CRR rescue approach remains available as a follow-on if deeper training alone is
insufficient.

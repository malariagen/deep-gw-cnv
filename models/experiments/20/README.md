# Experiment 20 — Curriculum-weighted VAE fine-tuning (cnv ramp 0.05→0.40)

**Status:** Complete 2026-04-19

**Full retraining required: VAE fine-tuned from exp 18 checkpoint with curriculum-weighted sampler.**

## Hypothesis

Exp 19 left GCH1 at FNR=0.07 (139 FNs, CRR p50=1.24, p90=1.29). The original exp 20
proposal was to lower `cnv_crr_amp_threshold` from 1.30 to 1.25 to recover FNs just below
the threshold. But feedback correctly pointed out that the bottleneck is not the threshold —
it is the VAE reconstruction quality. If the VAE has not learned CNV profiles well enough,
it reconstructs GCH1-amplified samples with near-normal accuracy, producing weak CRR signal
regardless of where the threshold sits.

The fix: retrain with a curriculum-weighted sampler. The sampler starts with a very low CNV
weight (ratio=0.05) so the VAE first consolidates its representation of normal read-count
profiles, then linearly ramps the CNV weight to 0.40 over 40 epochs. This exceeds the static
0.25 ratio used in experiments 18/19, giving the VAE substantially more exposure to CNV
profiles in the later phase of training. Better CNV representation → higher CRR for real
amplifications → GCH1 FNs cross the 1.30 threshold without any threshold change.

## Changes from experiment 19

| Component                         | Exp 19                   | Exp 20                            | Rationale                                              |
|-----------------------------------|--------------------------|-----------------------------------|--------------------------------------------------------|
| Sampler                           | fixed ratio=0.25         | curriculum 0.05→0.40 / 40 epochs  | More CNV exposure while protecting normal baseline     |
| `cnv_crr_amp_threshold`           | 1.30                     | 1.30 (unchanged)                  | Let VAE do the work; threshold stays at baseline       |
| VAE checkpoint (start)            | exp 18 symlink           | exp 18 (fine-tune from)           | Warm start from best existing checkpoint               |
| HMM segments                      | exp 18 symlink           | re-fitted from new VAE            | New VAE → new reconstructions → new segments           |
| Everything else                   | unchanged                | —                                 | Same architecture, sin loss, HMM params, CNV caller    |

## Expected outcome

- **GCH1 FNR**: improved VAE reconstructions should push CRR above 1.30 for samples that
  were at p50=1.24 in exp 19. Predicted FNR ≈ 0.03–0.05 if curriculum effectively boosts
  signal.
- **GCH1 PPV**: normal GCH1 p90=1.14 leaves a healthy gap below 1.30; PPV should stay
  ≥ 0.88.
- **PM2_PM3**: already at FNR=0.01; curriculum should not regress it. The fractional
  coverage caller + gate at 1.35 remains unchanged.
- **CRT / MDR1**: near-perfect and stable; not expected to change.

If GCH1 FNR barely improves, the signal weakness may be structural (population-specific
training distributions, particularly AF-E pre-2010 samples), and the next step would be
population-stratified sampling or per-population thresholds.

## Run time estimate

~2–3 hours (VAE training + HMM + CNV calling + evaluation on Mac mini).

## What we could do instead

1. **Population-specific thresholds** — AF-E (FNR=0.32) drives most of the GCH1 loss.
   A gene×population gate would concentrate the recovery without retraining.
2. **Lower threshold to 1.25** — the original exp 20 proposal; a quick test if retraining
   does not help. Only effective if signal is actually in [1.25, 1.30].
3. **self_transition=0.75** — still untested. Longer HMM segments may aggregate weak
   GCH1 signal better across the ~2 bins.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 0.86 | 0.73 |
| GCH1    | 0.05 | 0.71 | 0.79 |
| MDR1    | 0.01 | 1.00 | 0.93 |
| PM2_PM3 | 0.01 | 0.81 | 0.30 |

**GCH1**: FNR improved from 0.07 to 0.05 (98 FNs remain, CRR p50=1.25, p90=1.29 — barely
changed from exp 19's p50=1.24). The curriculum marginally boosted FN CRR but could not push
them above the 1.30 gate. More importantly, PPV collapsed from 0.91 to 0.71: 759 FPs with
CRR p50=1.35. The final ratio of 0.40 over-exposed the VAE to CNVs, inflating CRR for normal
samples past the gate. AF-E remains the worst population (FNR=0.14).

**PM2_PM3**: FNR unchanged at 0.01, but PPV regressed from 0.93 to 0.81 (22 FPs) — the same
curriculum inflation spilled into PM2_PM3. The prediction that PM2_PM3 would be unaffected
was wrong.

**CRT / MDR1**: Stable. MDR1 PPV improved to 1.00; CRT near-perfect.

**Divergence from predictions**: Both FNR and the impact on PPV diverged. FNR improved only
marginally (0.07→0.05 vs predicted 0.03–0.05), while PPV regressed severely (predicted
≥0.88). The curriculum at 0.40 inflated CRR globally rather than selectively boosting
GCH1-amplified signal.

**Next**: Experiment 21 — curriculum-weighted VAE fine-tuning with lower final ratio (0.05→0.30)
to reduce CRR inflation while retaining moderate CNV exposure improvement.

---

## Proposal history

### Original proposal (2026-04-19)

Lower `cnv_crr_amp_threshold` from 1.30 to 1.25 (GCH1 fallback gate). All 139 GCH1 FNs
had CRR p50=1.24, p90=1.29 — real signal just below the threshold. Predicted FNR drop
from 0.07 to 0.03–0.04. Reused exp 19 checkpoint and artifacts via symlinks; CNV calling
+ evaluation only (~30 min).

### Feedback received

Suggested exploring a curriculum WeightedRandomSampler that starts with a low CNV weight
(few CNVs seen early) and ramps up to 0.25–0.5, giving the VAE better CNV profile
exposure. Key point: if the raw CRR signal is low, the right fix is to boost the VAE,
not to adjust the HMM or call thresholds.

### What changed and why

The threshold tweak treats the symptom (CRR below the gate) rather than the cause (VAE
not reconstructing GCH1 amplifications strongly enough). Curriculum sampling addresses the
root cause. The approach is changed from a no-op symlink run to a full VAE retraining,
and the CNV threshold is reverted to the 1.30 baseline to test the VAE improvement in
isolation. The curriculum ends at ratio=0.40 (vs 0.25 static) to provide meaningfully
higher CNV exposure in the later training phase.

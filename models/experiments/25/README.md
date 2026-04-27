# Experiment 25 — Increase HMM self_transition 0.80→0.90 to reduce FPs from HMM jumpiness

**Status:** Complete 2026-04-22

**No retraining.** Reuses exp 21 checkpoint, latents, reconstructions, and sample_ids via
symlinks. Only the HMM, CNV caller, and evaluation are re-run.

## Hypothesis

Exp 24 produced p90=1167 within-chrom CN transitions per sample — an extreme level of
jumpiness for a haploid genome. This is most likely driving the FP surge seen across all
genes: CRT PPV 0.14 (354 FPs), MDR1 PPV 0.56 (552 FPs), PM2_PM3 PPV 0.41 (1227 FPs).
All three genes have FP CRR distributions starting at p10≈1.21, just above the 1.20 gate,
consistent with brief micro-calls in noisy regions that a stickier HMM would merge back
into CN1.

Raising hmm_self_transition from 0.80 to 0.90 increases the penalty for state switching,
forcing the Viterbi path to stay in a state longer before transitioning. The HMM
transitions count should fall substantially, cutting the spurious short-duration calls
without changing the threshold gate.

## Changes from experiment 24

| Component              | Exp 24 | Exp 25 | Rationale                                      |
|------------------------|--------|--------|------------------------------------------------|
| `hmm_self_transition`  | 0.80   | 0.90   | Reduce HMM jumpiness (p90=1167 transitions)    |
| Segments               | symlinked from exp 21 | re-run | Self_transition changed |
| Everything else        | unchanged | —   | Same 08 caller, 1.20 threshold, exp 21 model   |

## Expected outcome

- **CRT PPV**: 0.14 → ~0.25–0.40. FPs at CRR 1.21–1.25 are borderline micro-calls; a
  stickier HMM should merge them into the background CN1 state.
- **MDR1 PPV**: 0.56 → ~0.65–0.75. Similar mechanism.
- **GCH1 PPV**: 0.73 → ~0.80–0.87. Modest gain.
- **PM2_PM3 PPV**: 0.41 → ~0.50–0.60. Largest absolute FP count (1227); expect most
  reduction here.
- **FNR**: small increase possible (stickier HMM enters CN≥2 states less readily). Expect
  GCH1 FNR ≤ 0.03, CRT/MDR1 FNR = 0.00 unchanged.
- **PM2_PM3 FNs**: the 43 FNs with CRR p50=0.00 are coverage failures, not HMM issues;
  expect them to persist.
- **HMM transitions p90**: should fall from 1167 to well below 500.

## Run time estimate

~60–120 minutes (HMM fitting across 53,973 samples + CNV calling + evaluation; no training).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 0.14 | 0.37 |
| GCH1    | 0.02 | 0.73 | 0.83 |
| MDR1    | 0.00 | 0.56 | 0.73 |
| PM2_PM3 | 0.05 | 0.41 | 0.60 |

**HMM transitions p90: 1024** (down from 1167 in exp 24 — a 12% reduction, not the expected large drop).

Prediction vs reality:
- **PPV: no improvement at all.** All PPV values are identical to exp 24. The FP CRR distributions still start at p10=1.21 and cluster at 1.21–1.25 for every gene. The hypothesis that a stickier HMM would merge borderline micro-calls back into CN1 was wrong — either the signal is sustained rather than brief (so stickiness doesn't help) or these FPs come from genuine systematic CRR elevation near 1.20 rather than HMM jumpiness.
- **HMM transitions: modest reduction.** 1167 → 1024 (12%). Expected "well below 500". The HMM is still extremely jumpy; self_transition 0.90 is not sufficient to change its behaviour materially.
- **FNR: matched predictions.** GCH1 FNR stayed at 0.02, CRT/MDR1 at 0.00.
- **PM2_PM3 FNs: confirmed coverage failure.** 45 FNs with CRR p50=0.00 across all percentiles, unchanged.

The HMM tuning approach is exhausted at 0.90. The root cause of poor PPV is the 1.20 threshold admitting borderline CRR 1.21–1.25 calls that are systematically FPs. The fix is to raise the threshold, not tune the HMM. See experiment 26.

## Proposal summary

Proposed 2026-04-22. Single parameter change: hmm_self_transition 0.80→0.90. Motivated
by p90=1167 within-chrom transitions in exp 24 and FP CRR distributions starting at
1.21 for CRT, MDR1, and PM2_PM3. The stickier HMM was expected to reduce micro-calls in
noisy regions and recover PPV — but had no effect. Next: experiment 26 raises the
threshold back to 1.25.

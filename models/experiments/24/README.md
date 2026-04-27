# Experiment 24 — Lower amp threshold to 1.20 + fix PM2_PM3 GFF lookup + fixed-center logistic sharpness metric

**Status:** Complete 2026-04-22

**No retraining.** Reuses exp 21 checkpoint, latents, reconstructions, segments via symlinks.
Only the CNV caller and evaluation are re-run.

## Hypothesis

Exp 23 left 73 GCH1 false negatives, all with CRR in [1.07, 1.24] — clearly above the true-
normal p90 (1.13) but below the amp threshold (1.25). Lowering the threshold to 1.20 should
recover ~50% of these FNs (those with CRR ≥ 1.20) while preserving a 7pp buffer above 90%
of true normals.

PM2_PM3 was entirely absent from exp 23's evaluation because the caller used a merged
coordinate range (292244–299101) that matched no single GFF gene. In the PlasmoDB-54 GFF,
PM2 (PMII, PF3D7_1408000) occupies 292244–295261 and PM3 (PMIII, PF3D7_1408100) occupies
296683–299101. The fix: list both as separate REFERENCE_GENES entries with call_id "PM2_PM3",
then aggregate per-sample CN and CRR via max before pivoting.

The logistic sharpness metric (new in v08 caller) addresses a diagnostic gap: CRR alone does
not distinguish between a clean CN=2 segment with sharp boundaries and a noisy call where
individual bin values straddle the threshold. By fitting a sigmoid — with its center anchored
at the expected flank/gene boundary — to sorted per-bin copy ratios within the local window, we
get a per-call quality score that is robust to normalisation artefacts. A large `logistic_slope`
indicates that the transition in the sorted distribution occurs at the right structural position
(reliable CN boundary); a sluggish slope indicates noisy, overlapping distributions.

The fixed center is computed from the locus bin counts before NaN removal, so it reflects the
structural layout rather than the per-sample noise pattern. This prevents a flat sample with a
few elevated bins from scoring high sharpness by placing the sigmoid inflection wherever the
handful of noisy peaks happen to fall.

## Changes from experiment 23

| Component               | Exp 23 | Exp 24 | Rationale                                                      |
|-------------------------|--------|--------|----------------------------------------------------------------|
| `cnv_crr_amp_threshold` | 1.25   | 1.20   | Recover GCH1 FNs (CRR 1.07–1.24); buffer above normal p90=1.13 |
| PM2_PM3 GFF lookup      | broken | fixed  | PMII + PMIII listed separately; aggregated via max             |
| CNV caller version      | `06`   | `08`   | Adds `logistic_slope` with fixed-center logistic               |
| Everything else         | unchanged | —   | Same model, HMM, architecture, symlinks from exp 21            |

## Expected outcome

- **GCH1 FNR**: 0.04 → ~0.02. Expect to recover FNs with CRR ≥ 1.20 (~50% of the 73 FNs).
- **GCH1 PPV**: slight drop from 0.84, likely to ~0.78–0.82 (more borderline calls).
- **CRT/MDR1 PPV**: may decline modestly; FP CRR for both genes starts at p10=1.26, so
  calls at 1.20–1.25 are new territory — expect a small increase in FPs.
- **PM2_PM3**: now fully evaluated for the first time in the genome-wide pipeline.
- **logistic_slope**: amplified GCH1 calls should have higher slopes than the newly recovered
  borderline ones (CRR 1.20–1.24), providing a natural quality tier within the call set.

## Run time estimate

~60–120 minutes (CNV calling across 5,285 genes × 53,973 samples + logistic fits for all
amplified calls + 4 reference genes × all samples for diagnostics; no training or HMM).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| GCH1    | 0.02 | 0.73 | 0.83 |
| CRT     | 0.00 | 0.14 | 0.37 |
| MDR1    | 0.00 | 0.56 | 0.73 |
| PM2_PM3 | 0.05 | 0.41 | 0.60 |

**Predictions matched:**
- GCH1 FNR dropped from 0.04 → 0.02 exactly as predicted; HMM diagnostics showed 43 remaining FNs with CRR p50=1.17, consistent with soft signal below the new gate.
- PM2_PM3 evaluated for the first time — MCC=0.60.

**Predictions that diverged:**
- GCH1 PPV 0.73 was below the predicted 0.78–0.82; more borderline calls at CRR 1.20–1.24 than expected.
- CRT PPV 0.14 (worse than exp 23's 0.25): threshold lowering to 1.20 added a large wave of FPs at CRR 1.21–1.25 for CRT. MDR1 PPV similarly fell to 0.56 from 0.71.
- PM2_PM3 FNs had CRR p50=0.00 (all percentiles): not a threshold issue — these samples appear to have no gene signal relative to the flank, consistent with coverage failure or ground truth discrepancy.
- HMM produced p90=1167 within-chrom transitions per sample — extreme jumpiness not observed in prior experiments; likely driver of the FP surge.

**→ Experiment 25:** Increase HMM self_transition 0.80→0.90 to reduce jumpy micro-calls and improve PPV.

## Proposal history

**Original proposal (2026-04-21):** Lower `cnv_crr_amp_threshold` from 1.25 to 1.20 to
recover GCH1 FNs with CRR in [1.07, 1.24], and fix the PM2_PM3 GFF coordinate lookup.
No other changes. Used `06_genome_cnv_caller`.

**Feedback received (2026-04-22, round 1):** Add a logistic regression sharpness metric to the
CNV caller. For each amplified call, take the per-bin copy ratios in a local window around
the gene, sort them, and fit a sigmoid. A sharp gradient indicates a reliable CNV boundary;
a sluggish gradient indicates noisy data. Output this as a new column alongside CRR.
(Naming of the caller file left to discretion.)

**What changed (round 1):** A new caller version (`07_genome_cnv_caller`) added `_logistic_sharpness()`
which fits a sigmoid to sorted, normalised bin copy ratios within [gene_start - flank_padding,
gene_end + flank_padding]. The slope parameter was written as `logistic_slope` to both output
files. The threshold change and PM2_PM3 fix from the original proposal were unchanged.

**Feedback received (2026-04-22, round 2):** The free-center sigmoid is vulnerable to a
normalisation artefact: a flat sample with a few noisy peaks will have its y-axis compressed
to [0, 1], making those peaks look like a clean CN boundary. To prevent this, the center
should be fixed at the expected fraction of bins occupied by the flanking region, so that
only a transition occurring at the right structural position earns a high slope. Separately:
the top and bottom plateaux of the sigmoid should span lengths consistent with the gene and
flank bin counts.

**What changed (round 2):** Caller upgraded from `07_genome_cnv_caller` → `08_genome_cnv_caller`.
The sigmoid center is now fixed at `n_flank_in_window / n_total_in_window` (computed from bin
counts before NaN removal). Only the slope is fitted. This anchors the metric to the structural
layout of the locus: a genuine CNV boundary produces a sharp transition at the expected rank
position; a noisy sample cannot fake a high slope by shifting the inflection to wherever its
few peaks happen to fall. The threshold change and PM2_PM3 fix remain unchanged.

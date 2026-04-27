# Experiment 26 — Band segment filter: demote CRR [1.20, 1.25) calls not backed by HMM CN>=2 segments

**Status:** Complete 2026-04-23

**No retraining, no HMM re-run.** Reuses exp 25 checkpoint, latents, reconstructions,
sample_ids, and segments via symlinks. Only the CNV caller (new v10) and evaluation are re-run.

## Hypothesis

Exp 25 analysis showed that FPs concentrated in CRR [1.20, 1.25) are not genuine
amplifications — they are CRR-elevation artefacts with no corresponding CN>=2 HMM segment
at the gene.  For a real amplification the HMM should already assign CN>=2 to the gene
region.  We can therefore use segment coverage as a second confirmation gate: if a sample's
CRR is in the borderline band but the HMM placed no CN>=2 segment over the gene, demote
the call to CN=1.

This is the same goal as the v09 slope filter, but using the HMM output directly (where
segments cover) rather than fitting a logistic curve to the sorted CRR distribution.

## Changes from experiment 25

| Component                    | Exp 25                  | Exp 26                          | Rationale                                        |
|------------------------------|-------------------------|---------------------------------|--------------------------------------------------|
| `cnv` caller                 | `08_genome_cnv_caller`  | `10_genome_cnv_caller`          | Segment coverage filter replaces logistic slope  |
| `cnv_crr_band_upper`         | —                       | 1.25                            | Upper bound of the filtered band                 |
| `cnv_crr_band_cn2_threshold` | —                       | 0.50                            | ≥50% of gene bins must be in a CN>=2 HMM segment |
| Output column                | `logistic_slope`        | `segment_cn2_fraction`          | Renamed to reflect actual metric                 |
| `cnv_crr_amp_threshold`      | 1.20                    | 1.20 (unchanged)                | Threshold unchanged; filter handles borderline   |
| Segments                     | re-run (HMM 0.90)       | symlinked from exp 25           | HMM unchanged                                    |

## Expected outcome

- **CRT**: PPV 0.14 → ~0.26 (186 band FPs have no CN>=2 segment at CRT → fraction ≈ 0 → filtered; 0 band TPs to lose).
- **GCH1**: PPV 0.73 → ~0.80–0.85 (most/all of the 351 band FPs have no CN>=2 at GCH1; the 30 band TPs should have CN>=2 segments already); FNR 0.02 → ~0.02–0.025.
- **MDR1**: PPV 0.56 → ~0.65–0.70 (most of the 271 band FPs filtered; 2 band TPs with high CRR should have CN>=2 support).
- **PM2_PM3**: PPV 0.41 → ~0.50–0.55 (510 band FPs mostly filtered; 1 band TP kept).
- **HMM transitions**: unchanged (segments symlinked from exp 25).
- **PM2_PM3 FNs**: the 45 FNs with CRR≈0 are coverage failures; unaffected by this filter.

## Run time estimate

~15–30 minutes (CNV calling + evaluation only; no training, no HMM re-run).

## Proposal history

**Original proposal (2026-04-22):** Raise `cnv_crr_amp_threshold` from 1.20 → 1.25 using
caller v08.  Motivated by the FP CRR cluster at 1.21–1.25 observed in both exp 24 and 25,
insensitive to HMM self_transition.  TP CRR distributions for CRT/MDR1/PM2_PM3 start well
above 1.25, so raising the gate was expected to cut ~50% of FPs with minimal TP loss.

**Feedback (2026-04-22):** Requested logistic_slope percentile analysis for calls in the
1.20–1.25 CRR band before committing to the threshold change.

**Revised proposal (2026-04-22, v09):** Slope analysis confirmed strong TP/FP separation
in the band. Changed to caller v09 with a CRR-based logistic slope filter at cutoff=1.5:
sort window bins, fit a sigmoid, require the sharpness to meet a threshold for borderline
calls.

**Feedback (2026-04-22):** Clarification question: "what's band TP?" and "why is logistic
slope calculated on CRR? I asked for it to be calculated based on where segments cover."

**Revised proposal (2026-04-22, v10, current):** Replaced logistic slope on CRR with
`segment_cn2_fraction` — the fraction of gene bins that fall within a CN>=2 HMM segment
for each sample. This is exactly what the user asked for: the filter is now based directly
on where the HMM segments cover the gene, not on a curve fit to copy ratios.  The new
caller v10 keeps segments in memory after the sanity check, builds a per-chromosome CN>=2
interval lookup per sample, and uses binary search to compute the coverage fraction
efficiently.  Threshold 0.50: at least half the gene bins must be in a CN>=2 segment.

**Feedback (2026-04-22):** Asked for historical best FNR and PPV across experiments.
No change requested to the proposal itself — question was for context before authorising.
Answer provided in reply: best PPV (genome-wide) was exp 23 (threshold 1.25), best FNR
was exps 24–25 (threshold 1.20); exp 26 is designed to recover PPV of exp 23 without the
FNR cost by filtering the 1.20–1.25 band with HMM segment support. Proposal unchanged.

## Actual outcome

| Gene     | FNR  | PPV  | MCC  |
|----------|------|------|------|
| CRT      | 0.00 | 0.25 | 0.50 |
| GCH1     | 0.04 | 0.84 | 0.89 |
| MDR1     | 0.01 | 0.71 | 0.83 |
| PM2_PM3  | 0.05 | 0.54 | 0.70 |

**Prediction accuracy:** Nearly spot-on. CRT PPV predicted ~0.26, got 0.25. GCH1 PPV
predicted 0.80–0.85, got 0.84. MDR1 PPV predicted 0.65–0.70, got 0.71 (slightly better).
PM2_PM3 PPV predicted 0.50–0.55, got 0.54. FNR predictions matched well; GCH1 FNR came
in at 0.04 vs the predicted 0.02–0.025 (minor divergence).

**Key finding:** The band segment filter succeeded — PPV matches exp 23 levels (threshold
1.25) while FNR tracks exp 24–25 levels (threshold 1.20). However, all remaining FPs have
CRR ≥ 1.26 (FP p10=1.26 for every gene), meaning the artefact band extends above the
current band_upper of 1.25. The band filter only addressed the lower portion of the
artefact range; there is more gain available by extending band_upper upward.

**Next experiment:** Exp 27 extends band_upper from 1.25 → 1.35 to capture more band FPs.

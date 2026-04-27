# Experiment 27 — Extended band segment filter: widen artefact band from 1.25 → 1.35

**Status:** Complete 2026-04-23

**No retraining, no HMM re-run.** Reuses exp 25 checkpoint, latents, reconstructions,
sample_ids, and segments via symlinks. Only the CNV caller config and evaluation are re-run.

## Hypothesis

Exp 26 confirmed the band segment filter works: PPV recovered to exp 23 levels (threshold
1.25) while preserving the low FNR of exp 24–25 (threshold 1.20). But inspection of the
remaining FPs reveals a new key fact: **all remaining FPs have CRR ≥ 1.26** (FP p10=1.26
for every gene). The current band_upper of 1.25 means those FPs are entirely outside the
filter range — they pass through unchecked.

This tells us the artefact band extends above 1.25. The same artefact mechanism (systematic
CRR elevation without a real copy gain) continues into 1.25–1.35 and beyond. Extending
band_upper to 1.35 will subject those calls to the same CN>=2 HMM segment check that
successfully filtered the 1.20–1.25 artefacts in exp 26. Genuine amplifications in [1.25,
1.35) should retain HMM segment support and survive; artefacts without HMM backing will
be demoted to CN=1.

GCH1 is the most sensitive gene here: its TPs start at CRR p10=1.37 and p25=1.50, meaning
very few TPs fall in [1.25, 1.35). Those that do should have CN>=2 segment support (being
real amplifications). For CRT, MDR1, and PM2_PM3, TPs start at CRR p10 ≥ 1.67 — well
above 1.35 — so there is negligible TP risk.

## Changes from experiment 26

| Component             | Exp 26 | Exp 27 | Rationale                                             |
|-----------------------|--------|--------|-------------------------------------------------------|
| `cnv_crr_band_upper`  | 1.25   | 1.35   | Extend filtered band to cover more artefact-band FPs  |
| Everything else       | —      | —      | Unchanged                                             |

## Expected outcome

- **GCH1**: PPV 0.84 → ~0.88–0.92. GCH1 FPs have p75=1.33, so ~75% lie in [1.25, 1.35)
  and will be tested for HMM support; most are artefacts and will be demoted. GCH1 FNR
  0.04 → ~0.04–0.05 (slight risk from GCH1 TPs at CRR 1.37 that lack strong HMM segment
  coverage, but these should be HMM-backed real amplifications).
- **MDR1**: PPV 0.71 → ~0.79–0.83. MDR1 FPs have p75=1.34, so ~75% are in the new band.
  TPs start at p10=1.67 — no TP risk. FNR 0.01 → ~0.01.
- **PM2_PM3**: PPV 0.54 → ~0.61–0.66. PM2_PM3 FPs have p75=1.40, so ~50–60% lie below
  1.35 and will be filtered. FNR 0.05 → ~0.05 (PM2_PM3 FNs have CRR≈0, unaffected).
- **CRT**: PPV 0.25 → ~0.30–0.36. CRT FPs have p25=1.28 and p75=1.45; a meaningful
  fraction fall below 1.35. FNR 0.00 → 0.00.

## Run time estimate

~15–30 minutes (CNV calling + evaluation only; no training, no HMM re-run).

## Actual outcome

| Gene     | FNR  | PPV  | MCC  |
|----------|------|------|------|
| CRT      | 0.00 | 0.42 | 0.65 |
| GCH1     | 0.11 | 0.96 | 0.91 |
| MDR1     | 0.01 | 0.91 | 0.95 |
| PM2_PM3  | 0.06 | 0.73 | 0.82 |

**Prediction accuracy:** PPV gains exceeded every prediction. CRT: predicted 0.30–0.36,
got 0.42. GCH1: predicted 0.88–0.92, got 0.96. MDR1: predicted 0.79–0.83, got 0.91.
PM2_PM3: predicted 0.61–0.66, got 0.73. The band extension to 1.35 was more effective at
removing artefact FPs than expected.

**Key divergence — GCH1 FNR:** Predicted 0.04–0.05, got 0.11. The 223 GCH1 FNs have
CRR p50=1.28, p75=1.32, p90=1.34 — all concentrated in the new [1.25, 1.35) extended zone.
These are genuine amplifications at marginal CRR that lack sufficient CN>=2 HMM segment
coverage (below the 0.50 threshold). The single global threshold of 0.50 is too strict for
the extended zone where signal is weaker; the core [1.20, 1.25) zone likely has a cleaner
TP/FP separation. The remaining FPs (n=78) have CRR p10=1.34, meaning almost all sit above
band_upper=1.35 and bypassed the filter entirely.

**Next experiment:** Exp 28 introduces a two-tier band filter — preserving the 0.50 threshold
for the core [1.20, 1.25) zone while applying a more lenient 0.20 threshold in the extended
[1.25, 1.35) zone, to recover the GCH1 FNR regression without sacrificing the PPV gains.

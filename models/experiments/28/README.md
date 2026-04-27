# Experiment 28 — Two-tier band segment filter: core [1.20, 1.25) @ 0.50, extended [1.25, 1.35) @ 0.20

**Status:** Proposed 2026-04-23

**No retraining, no HMM re-run.** Reuses exp 25 checkpoint, latents, reconstructions,
sample_ids, and segments via symlinks. Only the CNV caller (new v11) and evaluation are re-run.

## Hypothesis

Exp 27 achieved outstanding PPV gains by extending the band filter to CRR 1.35, but at the
cost of GCH1 FNR rising from 0.04 to 0.11. The 223 GCH1 FNs have CRR p50=1.28 — all
concentrated in the newly filtered [1.25, 1.35) zone. These are genuine amplifications at
marginal CRR whose HMM segment support does not meet the global 0.50 threshold.

The root cause: at CRR 1.25–1.35, VAE reconstruction captures most of the signal and the
residual excess is small. The HMM at self_transition=0.90 may place the gene in a CN=2
segment only partially — enough to confirm a real copy gain exists, but not enough to cover
50% of gene bins. The core band [1.20, 1.25) has a cleaner separation (artefacts have
fraction≈0, TPs have fraction close to 1.0) so the 0.50 threshold works well there.

Caller v11 introduces a split-band design: the core zone keeps the proven 0.50 threshold,
while the extended zone [1.25, 1.35) uses a more lenient 0.20 threshold. This recovers TPs
with partial HMM coverage while still blocking pure artefacts (which should have fraction≈0
regardless of which CRR sub-range they fall in).

## Changes from experiment 27

| Component                      | Exp 27 | Exp 28 | Rationale                                                  |
|--------------------------------|--------|--------|------------------------------------------------------------|
| `cnv` caller                   | v10    | v11    | Two-tier band filter support                               |
| `cnv_crr_band_core_upper`      | —      | 1.25   | New: boundary between strict core and lenient extended zones |
| `cnv_crr_band_ext_cn2_threshold` | —    | 0.20   | New: lenient threshold for [1.25, 1.35) extended zone      |
| `cnv_crr_band_cn2_threshold`   | 0.50   | 0.50   | Unchanged: strict threshold for [1.20, 1.25) core zone     |
| `cnv_crr_band_upper`           | 1.35   | 1.35   | Unchanged: outer boundary                                  |

## Expected outcome

- **GCH1**: FNR 0.11 → ~0.05–0.07. TPs with partial HMM coverage (frac 0.20–0.49) in the
  extended zone are now recovered. PPV 0.96 → ~0.91–0.95 (slight regression as some artefacts
  with modest HMM overlap pass through the extended zone; but remaining artefacts should be
  far fewer than the TPs recovered).
- **MDR1**: FNR 0.01 → ~0.01 (unchanged). MDR1 FNs have CRR p50=1.23 — in the core zone
  where the threshold is unchanged. MDR1 TPs all have CRR ≥ 1.67. No impact expected.
- **PM2_PM3**: FNR 0.06 → ~0.06 (unchanged). FNs have CRR p50=0.00 (coverage failures
  not related to the band filter).
- **CRT**: FNR 0.00 → 0.00, PPV 0.42 → ~0.42 (unchanged). All remaining FPs have
  CRR ≥ 1.36 — above the extended zone's upper bound; unaffected by this change.

## Run time estimate

~15–30 minutes (CNV calling + evaluation only; no training, no HMM re-run).

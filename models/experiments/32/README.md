# Experiment 32 — Lower HMM self_transition 0.96 → 0.90 to recover gene-scale FNs

**Status:** Complete 2026-04-29

**No retraining, no re-inference.** Reuses exp 30 checkpoint, latents, reconstructions,
and sample_ids via symlinks. HMM is re-run at the new self_transition value; segments
cannot be reused from exp 31.

## Hypothesis

Exp 31 switched to HMM-segment-primary calling (v14 caller, conf=0.25) and achieved a
PPV improvement, but FNR regressed severely: GCH1 0.03→0.49, MDR1 0.00→0.24,
PM2_PM3 0.02→0.46, CRT 0.02→0.12.

The near-miss diagnostics are unambiguous: 87–95% of FNs have `segment_cn2_fraction=0`
at confidence 0.25. FN CRR p50 values are 1.87–2.73 — these are genuine amplifications.
The HMM is not producing any CN>=2 coverage at these gene loci, so no threshold tuning
can recover them. This is an HMM segmentation problem, not a calling-threshold problem.

The evaluation guidance says exactly this: "FN p50 >> 1.0 → signal present, HMM
discarding it. Fix: lower self_transition."

At `self_transition=0.96`, the expected segment duration is ~25 bins × 400bp = ~10kb.
Gene-scale amplifications may span only 7–15 bins (3–6kb) at 400bp resolution. If the
elevated-CRR signal doesn't persist long enough to match the HMM's prior for a long
segment, the HMM stays in CN=1. Lowering to `self_transition=0.90` reduces expected
segment duration to ~10 bins × 400bp = ~4kb — appropriate for gene-scale events.

Note: exp 29 used st=0.90 at 400bp and showed p50=270 transitions ("over-segmentation"),
which motivated the 0.96 recalibration in exp 30. The concern then was fragmented segments
failing the minimum-coverage filter. With HMM-primary calling and coverage_fraction=0.50,
shorter but denser segments should still pass the majority-of-gene-body criterion. The
near-miss diagnostics will confirm whether frac=0 FNs are recovered.

## Changes from experiment 31

| Component              | Exp 31                     | Exp 32                     | Rationale                                              |
|------------------------|----------------------------|----------------------------|--------------------------------------------------------|
| `hmm_self_transition`  | 0.96                       | **0.90**                   | Shorter expected segments recover gene-scale FNs       |
| Segments               | symlinked from exp 30      | **re-run**                 | New st value requires new HMM run                      |
| `out_dir`              | `..._st0.96_...`           | **`..._st0.90_...`**       | Reflects the self_transition value                     |
| Everything else        | same                       | same                       |                                                        |

## Expected outcome

- **FNR**: large reduction expected. Most FNs have frac=0 because the HMM stays in
  CN=1 at st=0.96. At st=0.90, shorter segments are allowed; genuine amplifications
  with 7–15 elevated bins should produce CN>=2 segments that pass the 50% threshold.
  GCH1 is the most uncertain: only ~4–5 bins at 400bp resolution, so even st=0.90
  may be marginal. MDR1 and PM2_PM3 have more bins and should recover strongly.

- **PPV**: the main risk is that lowering st produces fragmented CN>=2 segments on
  normal samples too, increasing FPs. FP CRR p50 values from exp 31 are 1.29–1.40
  (near normal), suggesting FP samples have weak or absent HMM signal. If so, the
  50% coverage_fraction filter will block them even with shorter segments. We should
  monitor FP `segment_cn2_fraction` distributions in the near-miss diagnostics.

- **Transitions per sample**: p50 should increase from ~215 (st=0.96) toward ~270
  (st=0.90, as seen in exp 29). This is expected and acceptable.

- **MCC**: net improvement expected if FNR recovery dominates any PPV regression.

## Run time estimate

~2–3 hours (HMM re-run on 50K samples ~2h; CNV calling + evaluation ~30 min).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.11 | 0.14 | 0.35 |
| GCH1    | 0.45 | 0.70 | 0.58 |
| MDR1    | 0.21 | 0.53 | 0.63 |
| PM2_PM3 | 0.39 | 0.63 | 0.61 |

Lowering st from 0.96 to 0.90 produced only marginal FNR improvement (3–7pp vs exp 31),
far short of the predicted "substantial recovery." PPV also regressed slightly (as
predicted from increased fragmentation). The transitions p50 was 270, exactly as
expected from exp 29 — the HMM IS more active, but not at the gene loci.

The near-miss diagnostics are the key finding: 93% of GCH1 FNs, 62% of MDR1 FNs, and
83% of PM2_PM3 FNs still have `segment_cn2_fraction=0` at st=0.90. FN CRR p50 values
remain high (GCH1=1.84, MDR1=2.66, PM2_PM3=2.53) — the signal is real, but the HMM
cannot segment it. The root cause is that gene-specific amplifications spanning only 4–5
bins are overwhelmed in the per-sample HMM fit by the rest of the chromosome; the CN=2
Gaussian state fits to stronger signals elsewhere, leaving the gene-level 1.5–2.0
copy-ratio elevation classified as CN=1. Lowering st addresses transition probability
but not the fitting problem. Further st reductions are not the right fix.

→ Experiment 33: introduce a CRR rescue fallback in the caller (v15) — call CN=2 if
HMM misses but CRR ≥ 1.60. No HMM re-run needed; reuses exp 32 segments.

# Experiment 31 — HMM-segment-primary calling, confidence 0.25 (v14 caller)

**Status:** Complete 2026-04-28

**No retraining, no HMM re-run.** Reuses exp 30 checkpoint, latents, reconstructions,
sample_ids, and segments.parquet via symlinks. Only CNV calling and evaluation are re-run.

## Hypothesis

Exp 30 achieved near-zero FNR across all genes but PPV collapsed: CRT=0.07, MDR1=0.52,
PM2_PM3=0.50. The root cause was misidentified in the first two proposals: the FP problem
is not a CRR threshold problem — it is an architectural problem. CRR was designed as a
diagnostic signal indicator, not as the primary calling criterion. Using it as the
primary gate, with a supplementary 20% CN>=2 HMM coverage filter in the "band" zone, is
scientifically indefensible: 20% coverage means only 1/5 of the gene body is in CN>=2,
which is not a principled claim that the gene is amplified.

v13 discards CRR as a gate entirely. The primary call criterion is:

    >= 50% of the gene's bins are covered by high-confidence CN>=2 HMM segments.

This is scientifically defensible: the majority of the gene body must be in an elevated
copy-number state. The chromosome-level sanity check (>= 55% CN=1 background on the
chromosome) is retained — it validates that the HMM has correctly established the CN=1
baseline before making any gene-level call.

CRR is still computed (using the same low-coverage filter from v12) and stored in the
output as a diagnostic column; it does not influence the call decision.

## Changes from experiment 30

| Component                    | Exp 30                     | Exp 31                         | Rationale                                         |
|------------------------------|----------------------------|--------------------------------|---------------------------------------------------|
| `cnv`                        | `12_genome_cnv_caller`     | **`14_genome_cnv_caller`**     | HMM-primary + lower confidence + CN>=3 explicit   |
| `evaluation`                 | `03_pf9_evaluation`        | **`04_pf9_evaluation`**        | Adds near-miss / blindspot diagnostics            |
| `cnv_min_confidence`         | 0.50                       | **0.25**                       | Include lower-confidence segments to recover FNs  |
| `cnv_min_gene_coverage_fraction` | 0.50 (unused in v12)   | **0.50 (primary call criterion)** | ≥50% of gene body in CN>=2 is defensible       |
| CN>=3 handling               | implicit via `cn>=2`       | **explicit in v14 naming/docs** | CN=2→CN=3 transitions within gene body correct   |
| CRR band parameters          | `band_upper`, `band_cn2_threshold`, etc. | **removed** | No longer used for calling           |
| `out_dir`                    | `..._cnvcaller_fix`        | **`..._hmm_primary_conf0.25`** | Reflects approach and confidence level            |
| Segments                     | computed                   | symlinked from exp 30          | HMM unchanged; no re-run needed                  |
| Everything else              | same                       | same                           |                                                   |

## Expected outcome

- **PPV**: substantial recovery across all genes. The FPs from exp 30 (CRR 1.35–1.50, passing
  unconditionally) almost certainly do not have 50% CN>=2 gene-body coverage — noise-driven
  CRR elevations are driven by a handful of bins and produce short, fragmented HMM events.
  CRT PPV: 0.07 → likely 0.40–0.70. MDR1 PPV: 0.52 → likely 0.70–0.90.

- **FNR**: lowering confidence to 0.25 should recover FNs caused by low-confidence HMM
  segments at genuine amplifications. The primary risk is that low-confidence CN>=2
  segments are noise-driven; if they are, PPV will degrade. GCH1 is most at risk from
  noise: ~4–5 bins at 400bp, so borderline segments could swing the call either way.
  GCH1 FNR: uncertain, could improve (0.03 → ~0.01) or stay flat if GCH1 FNs are already
  well-covered at higher confidence.

- **Near-miss diagnostics** (new in v04 evaluator): will show FN segment_cn2_fraction
  distributions and recovery counts, making the next iteration decision concrete.

- **MCC**: net improvement expected if PPV recovery outweighs any FNR regression.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.12 | 0.18 | 0.39 |
| GCH1    | 0.49 | 0.75 | 0.58 |
| MDR1    | 0.24 | 0.58 | 0.65 |
| PM2_PM3 | 0.46 | 0.65 | 0.58 |

**Where predictions matched:** PPV improved across all genes (CRT 0.07→0.18, GCH1 0.63→0.75,
MDR1 0.52→0.58, PM2_PM3 0.50→0.65), confirming that the HMM-primary criterion is more
precise than the old CRR-gated band approach.

**Where predictions diverged:** FNR regressed dramatically on all genes, especially GCH1
(0.03→0.49) and PM2_PM3 (0.02→0.46). The expectation was that lowering confidence to 0.25
would recover FNs — it did not. Near-miss diagnostics show the cause: 87–95% of FNs have
`segment_cn2_fraction=0`, meaning the HMM assigned zero CN>=2 coverage to these samples at
their gene loci at any confidence level. FN CRR p50 values are 1.87–2.73, confirming the
signal is present in raw coverage. The HMM at `self_transition=0.96` (expected segment
duration ~10kb) is simply not segmenting gene-scale amplifications. Threshold tuning alone
cannot fix a frac=0 FN — the HMM must produce CN>=2 segments first.

**Next:** Experiment 32 — lower `hmm_self_transition` from 0.96 → 0.90 to reduce expected
segment duration from ~10kb to ~4kb, making the HMM responsive to gene-scale events.

## Run time estimate

~30 minutes (CNV calling + evaluation only; segments already exist from exp 30).

## Proposal history

**First proposal (2026-04-28):** Raise `cnv_crr_band_upper` from 1.35 → 1.55 to gate
the FP-dense 1.35–1.55 CRR range through the 20% HMM CN>=2 check. No HMM re-run; reuse
exp 30 segments via symlink.

**Feedback after first proposal (2026-04-28):** Asked whether GCH1 is still being treated
differently from other genes and whether 400bp resolution makes GCH1 no longer CRR-dependent.

**Second proposal (2026-04-28):** Parameter change unchanged (band_upper 1.35 → 1.55).
Clarified that GCH1 receives no special treatment in v12 (or any caller from v06 onwards),
and that the dead `cnv_crr_min_bins_fallback` config key was removed to avoid confusion.

**Feedback after second proposal (2026-04-28):** CRR was never intended to be the primary
calling criterion; it was always only a diagnostic signal indicator. The 20% CN>=2 coverage
threshold in the extended band zone is scientifically indefensible.

**Third proposal (2026-04-28):** Dropped the entire CRR-gated band filter. Introduced
HMM-segment-primary caller (v13): gene called amplified if >=50% of bins covered by
high-confidence (confidence>=0.50) CN>=2 HMM segments. CRR retained as diagnostic only.

**Feedback after third proposal (2026-04-28):** Lower confidence threshold to 0.25 (from
0.50) to include more segments. Also ensure CN=2→CN=3 transitions within the gene body are
handled well. Additionally, design new evaluation diagnostics for blindspots and near-misses;
the evaluation suite should be treated as a living document — add new diagnostics when they
expose specific failure modes, deprecate metrics that are no longer relevant, and always
require scientific defensibility.

**What changed in this (fourth) proposal:**
- Caller updated to v14: same HMM-primary logic but with explicit CN>=3 handling
  (renamed internals, added documentation confirming CN=2→CN=3 mid-gene transitions are
  handled by unioning all elevated-CN intervals).
- `cnv_min_confidence` lowered from 0.50 → 0.25 in config (no code change needed).
- Evaluation updated to v04: adds NEAR-MISS DIAGNOSTICS section showing, for each FN,
  the distribution of `segment_cn2_fraction` and explicit recovery counts at threshold
  0.35 and 0.25. Complete-miss FNs (fraction=0) are flagged separately.
- `out_dir` renamed to `..._hmm_primary_conf0.25` to reflect the confidence level.

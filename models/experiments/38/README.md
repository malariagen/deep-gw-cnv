# Experiment 38 — Upgrade to 05_pf9_evaluation using exp 37 artifacts

**Status:** Complete 2026-05-13

**No retrain, no HMM re-run, no threshold change.** Symlinks exp 37 checkpoint, reconstructions,
and segments. Single change: upgrade from `04_pf9_evaluation` to `05_pf9_evaluation`.

## Hypothesis

Experiment 37 achieved GCH1 FNR=0.09 (target ≤0.10 met). The evaluation metric (FNR) only covers
samples where both ground truth and model are callable — samples the model calls as -1 are silently
excluded before FNR is computed. If ground-truth CNV samples fail quality filters at a
disproportionately high rate, FNR understates total recall loss and we cannot see it.

The `05_pf9_evaluation` adds a **masked CNV diagnostic** to make this blindspot visible before
any threshold changes are made, establishing a clean baseline.

## Changes from experiment 37

| Parameter    | Exp 37                | Exp 38                | Rationale                                    |
|--------------|-----------------------|-----------------------|----------------------------------------------|
| `evaluation` | `04_pf9_evaluation`   | **`05_pf9_evaluation`** | Adds masked CNV rate and total FNR metrics  |

All other parameters (including `cnv_crr_rescue_threshold` = 1.55) and all artifacts unchanged.
Artifacts symlinked from exp 37.

## Expected outcome

- **FNR / PPV / MCC**: identical to exp 37 — no calling parameters changed.
- **masked_cnv_rate**: expected to be low (≈ overall model_miss of ~0.02) if -1 calls are unbiased.
  If substantially higher, the calling pipeline disproportionately suppresses CNV samples.
- **total_fnr**: expected to closely match FNR if masked_cnv_rate is low.
- **Masked CRR distribution**: provides signal on whether masked CNV samples are rescuable
  (high p50 = they have detectable signal; low p50 = they are genuinely low coverage).

This establishes the baseline needed to safely interpret any future threshold or HMM changes.

## Evaluation upgrade (05_pf9_evaluation)

New outputs:
- `masked_cnv_rate` — fraction of gt-positive samples called as -1 per gene
- `total_fnr` — (FN_in_eval + masked_cnv_n) / gt_pos_n; inclusive miss rate
- `MASKED CNV CRR` — CRR distribution for masked CNV samples

If `total_fnr ≈ FNR`, the -1 calls are unbiased. If `total_fnr >> FNR`, the calling pipeline
is disproportionately suppressing CNV samples — a blindspot that FNR alone cannot detect.

## Run time estimate

~15 minutes (CNV calling + evaluation only; all artifacts symlinked from exp 37).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 1.00 | 1.00 |
| GCH1    | 0.09 | 1.00 | 0.94 |
| MDR1    | 0.02 | 1.00 | 0.99 |
| PM2_PM3 | 0.00 | 1.00 | 1.00 |

**Predictions partially matched — but the key finding was unexpected.**

FNR/PPV/MCC were identical to exp 37 as expected (no calling parameters changed). However,
masked_cnv_rate was catastrophically high — 84–98% of ground-truth CNV samples were called -1:

| Gene    | gt_pos | masked_n | masked_rate | total_fnr |
|---------|--------|----------|-------------|-----------|
| CRT     | 58     | 57       | 0.98        | 0.98      |
| GCH1    | 2036   | 1815     | 0.89        | 0.90      |
| MDR1    | 717    | 603      | 0.84        | 0.84      |
| PM2_PM3 | 949    | 885      | 0.93        | 0.93      |

**The expected masked_cnv_rate was ~0.02 (same as overall call_rate).** Instead, 89–98% of
CNV-positive samples are suppressed before evaluation. The masked CRR is high (GCH1 p50=2.10,
all genes p50 >> 1.0), confirming the signal IS present — quality filters are blocking the calls.

**Root cause identified:** `cnv_max_transitions: 50` in the noise filter (v16 caller) is the culprit.
The HMM produces transitions p10=85, p50=1189 per sample — essentially ALL samples exceed 50
transitions. The noise filter fires at priority 0, before sanity check and before CRR rescue,
so even CNV samples with CRR=3.0 are called -1.

**Next experiment:** Exp 39 — disable the noise filter (`cnv_max_transitions: null`) to stop
silently masking CNV samples with legitimate high-CRR signal.

## Proposal history

**First proposal (2026-05-13):** Lower `cnv_crr_rescue_threshold` 1.55 → 1.35 using exp 37
artifacts. No other changes. Evaluation version remained at `04_pf9_evaluation`.

**Feedback 1:** Asked whether samples with ground-truth CNVs are being called as -1, and whether
the evaluation has a metric to catch this. The answer was no — v04 FNR only covers the evaluable
subset, so ground-truth CNV samples the model calls as -1 are silently excluded.

**Second proposal (2026-05-13):** Keep the threshold change (1.55 → 1.35) AND upgrade to
`05_pf9_evaluation`. Both changes bundled together.

**Feedback 2:** Only improve evals this experiment — no threshold change.

**Revised proposal:** Drop the threshold change entirely. Exp 38 only upgrades to `05_pf9_evaluation`
to establish the masked CNV baseline. The threshold change (1.55 → 1.35) can be evaluated as a
separate experiment once we understand the -1 suppression rate.

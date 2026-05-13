# Experiment 35 — Lower CRR rescue threshold (1.60 → 1.45) to recover GCH1 FNs

**Status:** Complete 2026-05-06

**No retraining, no new HMM run.** Symlinks exp 32 checkpoint, reconstructions, and
segments. Single config change: `cnv_crr_rescue_threshold` 1.60 → 1.45.

## Hypothesis

Exp 34 established that the noise filter + CRR rescue architecture works: PPV improved
dramatically to 0.89–0.99, call_rate ≈ 13% (consistent with the ~15% prediction), and
CRR rescue fired for 593 GCH1 TPs, 317 MDR1 TPs, and 254 PM2_PM3 TPs.

The remaining problem is GCH1 FNR=0.20 (147 FNs). Near-miss diagnostics show 146/147 have
`segment_cn2_fraction=0` — the HMM gave no CN>=2 coverage for these samples. But CRR BY
CALL OUTCOME shows these FNs have genuine CRR signal (p25=1.38, p50=1.48, p90=1.58) that
sits just below the 1.60 rescue threshold.

The threshold was originally calibrated against exp 32's noisy full-sample TN background
(TN CRR p90=1.50–1.53). After the noise filter, the callable subset has a much tighter TN
distribution (p90=1.12). This means 1.60 is now overly conservative: there is a large gap
between TN CRR p90=1.12 and the FN CRR p25=1.38. Lowering the threshold to 1.45 sits in
this gap and should recover the majority of GCH1 FNs.

## Changes from experiment 34

| Parameter                  | Exp 34  | Exp 35   | Rationale                                              |
|----------------------------|---------|----------|--------------------------------------------------------|
| `cnv_crr_rescue_threshold` | 1.60    | **1.45** | Callable TN p90=1.12; FN p25=1.38 → safe to lower     |

All other parameters and all artifacts unchanged.

## Expected outcome

- **GCH1 FNR**: should drop from 0.20. At threshold 1.45, expect to rescue FNs with
  CRR ≥ 1.45 ≈ FN p25 (roughly 75%+ of 147 FNs recovered). Target FNR < 0.10.
- **PPV**: should hold. In the callable subset, TN CRR p90=1.12 — no TN density between
  1.12 and 1.45. The existing 5 GCH1 FPs have CRR p10=1.28, so a small number of new FPs
  from the 1.28–1.44 CRR range is possible but count should remain small.
- **Other genes**: MDR1/CRT/PM2_PM3 FNR already near zero; effect there is minimal.
- **Missingness and call_rate**: unchanged (noise filter is unchanged).

## Run time estimate

~30 minutes (CNV calling + evaluation only; all artifacts symlinked from exp 32).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| GCH1    | 0.09 | 0.97 | 0.93 |
| MDR1    | 0.01 | 0.88 | 0.93 |
| CRT     | 0.04 | 0.55 | 0.72 |
| PM2_PM3 | 0.00 | 0.81 | 0.90 |

call_rate = 13% for all genes (noise filter unchanged).

**GCH1 FNR prediction matched.** Target was < 0.10; achieved 0.09. 66 FNs remain (down
from 147). All 66 have segment_cn2_fraction=0 and CRR p90=1.43 — just below the 1.45
threshold. GCH1 PPV held at 0.97 as expected.

**CRT PPV collapsed unexpectedly.** In exp 34 (threshold 1.60) CRT PPV was 0.92. Lowering
to 1.45 introduced 16 new CRT FPs with CRR 1.45–1.52 (p10=1.45, p90=1.52). CRT TPs
remain at CRR ≥ 1.66 (p10=1.66). The global threshold optimised for GCH1 is too low for
CRT — these two genes need separate thresholds.

**PM2_PM3 PPV also regressed.** PPV dropped 0.90 → 0.81 (+30 FPs with CRR 1.46–1.63).
Same cause: global threshold 1.45 is too low for PM2_PM3.

**Next experiment:** Exp 36 — per-gene CRR rescue thresholds (GCH1 1.35 / CRT 1.55 /
MDR1 1.45 / PM2_PM3 1.50), using new v17 caller.

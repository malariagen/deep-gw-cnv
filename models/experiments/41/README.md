# Experiment 41 — Raise HMM self_transition 0.96 → 0.99 to fix sanity-check failure

**Status:** Complete 2026-05-13

**No retrain, no inference re-run.** Reuses exp 37 checkpoint, latents, and reconstructions
via symlinks. HMM segmentation is re-run with the new self_transition. CNV calling and
evaluation are re-run.

## Hypothesis

Experiment 40 raised `hmm_self_transition` from 0.90 to 0.96, but transitions dropped
only from p50=1189 to p50=970 — an 18% reduction, far below the expected O(10–100). The
sanity check (≥55% chromosome CN=1) continues to mask 40–67% of GT CNV positives, all of
which have high CRR (p25=1.75–2.00, p50=1.76–2.44 across genes) confirming signal is present.

**Root cause refinement:** The 06_conv_vae reconstruction residuals are fundamentally noisier
than the VAE used in exp 30/31. At self_transition=0.96, the HMM state-transition prior
(1-p=0.04) is not tight enough to suppress the rapid oscillations driven by the new VAE's
emission variance. The relationship between self_transition and transitions is less linear
than expected — simply moving from 0.90 to 0.96 was insufficient.

**Proposed fix:** Raise `hmm_self_transition` from 0.96 to 0.99. At 0.99, (1-p)=0.01 vs
0.04 at 0.96 — 4× lower transition probability. Empirical extrapolation from the 0.90→0.96
observations suggests p50 transitions of ~240, which is in the range of exp 30 (p50=270)
where sanity checks passed. CRR rescue at 1.55 remains the safety net for amplifications
the stickier HMM fails to segment (masked CRR p25=1.75–2.00 in exp 40).

## Changes from experiment 40

| Parameter              | Exp 40 | Exp 41   | Rationale                                                                    |
|------------------------|--------|----------|------------------------------------------------------------------------------|
| `hmm_self_transition`  | 0.96   | **0.99** | 4× lower transition probability; targets p50~240 to restore sanity-check pass |
| HMM re-run             | yes    | **yes**  | Self-transition change requires new segments.parquet                         |

All other parameters unchanged. Noise filter remains disabled (`cnv_max_transitions: null`).

## Expected outcome

- **HMM transitions p50**: drops from 970 toward ~240 (extrapolated from 0.90→0.96 ratio).
- **Sanity check pass rate**: most samples pass, collapsing the masked CNV pool.
- **masked_cnv_rate**: drops from 40–67% toward a small residual.
- **CRR rescue role**: with a stickier HMM, more FNs will have segment_cn2_fraction=0 but
  CRR rescue (≥1.55) catches them — masked CRR p25=1.75–2.00 across all genes is well above
  the threshold.
- **total_fnr**: should drop toward the evaluable FNR (0.02–0.11 from exp 40).
- **PPV**: the evaluable pool grows as more samples clear sanity; PPV may shift slightly.
  TN CRR p90=1.07–1.24 in exp 40 — very few true negatives near the 1.55 rescue threshold.
- **CRT**: still challenging (58 GT positives, high FP rate, low PPV in all experiments).

## Run time estimate

~2–3 hours: HMM re-segmentation (~2h with n_jobs=-1 across 50k samples) +
CNV calling (~30 min) + evaluation (~5 min). No training or inference.

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.06 | 0.02 | 0.12 |
| GCH1    | 0.10 | 0.62 | 0.71 |
| MDR1    | 0.02 | 0.59 | 0.75 |
| PM2_PM3 | 0.04 | 0.39 | 0.59 |

Masked CNV rate: CRT=0.71, GCH1=0.47, MDR1=0.40, PM2_PM3=0.52.
Total FNR: CRT=0.72, GCH1=0.53, MDR1=0.42, PM2_PM3=0.54.
HMM transitions p50=769 (expected ~240; actual only 21% reduction from exp 40's 970).

**Where predictions matched:** Evaluable FNR remains excellent (0.02–0.10). Masked CRR
confirms signal is present (p10 ≥ 1.48 across all genes), and near-miss diagnostics
confirm nearly all FNs have segment_cn2_fraction=0 — HMM not segmenting them at all.

**Where predictions diverged:** The expected 4× reduction in (1-p) (0.04 → 0.01) did not
translate into a 4× reduction in transitions — only 21%, continuing the pattern from exp 40.
The masked_cnv_rate barely changed vs exp 40 (CRT actually worsened: 0.67 → 0.71). This
confirms the issue is not the transition prior but the emission noise: the 06_conv_vae
residuals have intrinsically high bin-to-bin variance that the self_transition prior cannot
suppress. Raising self_transition further has hit diminishing returns.

**Next:** Experiment 42 — smooth HMM emission input (rolling mean, 5-bin window) to reduce
bin-to-bin copy ratio variance before the HMM sees it.

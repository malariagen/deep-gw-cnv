# Experiment 42 — Retrain 06_conv_vae without sin loss: clean 400bp baseline

**Status:** Complete 2026-05-14

**Full retrain.** Sin loss removal changes the training objective, so no prior checkpoint
can be reused. HMM segmentation uses `03_gaussian_hmm` with no smoothing — this is the
clean 400bp baseline.

## Hypothesis

Experiment 41 confirmed that three successive raises of `hmm_self_transition` (0.90 → 0.96 → 0.99)
reduced HMM transitions p50 by only 35% in total (1189 → 970 → 769), and masked CNV rates
remained 40–71% across genes. The prior diagnosis was that 06_conv_vae residuals have high
per-bin variance, overwhelming the HMM's transition prior.

**Root cause:** The sin loss term (`sin_loss_max_weight=2.0`) introduced in exp 37 to
boost GCH1 reconstruction sensitivity is producing noisy, weird reconstructions as a side effect.
The sinusoidal regulariser adds an oscillatory component to the loss landscape that translates into
per-bin reconstruction artefacts — exactly the kind of noise that makes adjacent CN=1 bins look
like CN=2 to the HMM. Removing sin loss should yield cleaner residuals, making CN states more
separable and reducing spurious HMM transitions.

## Changes from experiment 41

| Parameter                 | Exp 41 | Exp 42    | Rationale                                                                  |
|---------------------------|--------|-----------|----------------------------------------------------------------------------|
| `sin_loss_max_weight`     | 2.0    | **0.0**   | Sin loss causing noisy reconstructions; removing it to clean VAE residuals |
| `sin_loss_warmup_epochs`  | 30     | **0**     | No longer relevant with max_weight=0.0                                     |
| VAE retrain               | no     | **yes**   | Training objective changed; exp 37 checkpoint no longer valid              |
| `hmm` module              | `03_gaussian_hmm` | `03_gaussian_hmm` | Unchanged — clean baseline, no rolling-mean smoothing      |
| `hmm_self_transition`     | 0.99   | 0.99      | Unchanged                                                                  |

All other parameters unchanged. Noise filter remains disabled (`cnv_max_transitions: null`).

## Expected outcome

- **VAE reconstructions**: smoother without sinusoidal artefacts; per-bin residual variance
  should drop substantially, making CN=1 and CN=2 distributions more separable.
- **HMM transitions p50**: should drop well below 769 — target <200. Cleaner emissions reduce
  the number of bins where the HMM sees a false state-change opportunity.
- **Sanity check pass rate**: most samples should pass (≥55% chromosome CN=1), collapsing the
  masked CNV pool.
- **masked_cnv_rate**: drops from 40–71% toward a small residual.
- **CRR rescue role**: for FNs where the HMM still fails to segment, CRR rescue (≥1.55) fires.
  Masked CRR p10 ≥ 1.48 across all genes in exp 41 — most masked CNV samples rescuable.
- **total_fnr**: should converge toward evaluable FNR (0.02–0.10 from exp 41).
- **GCH1 sensitivity risk**: sin loss was added in exp 37 to help GCH1 in AF-W/AF-E where
  frac=0 failures occurred. Without sin loss, GCH1 FNR may worsen in Africa; worth monitoring
  closely in by-population results.

## Run time estimate

~7–8 hours: VAE training (~5–6h for 100 epochs on 50k samples) + HMM re-segmentation (~2h)
+ CNV calling (~30 min) + evaluation (~5 min).

## Actual outcome

| Gene    | FNR  | PPV  | MCC  |
|---------|------|------|------|
| CRT     | 0.00 | 0.05 | 0.21 |
| GCH1    | 0.14 | 0.44 | 0.55 |
| MDR1    | 0.00 | 0.25 | 0.48 |
| PM2_PM3 | 0.00 | 0.17 | 0.37 |

Masked CNV rate: CRT=0.12, GCH1=0.17, MDR1=0.19, PM2_PM3=0.33.
Total FNR: CRT=0.12, GCH1=0.29, MDR1=0.19, PM2_PM3=0.33.
HMM transitions p50=306 (target was <200; actual 60% reduction from exp 41's 769).
Callability: 0.995.

**Where predictions matched:** Sin loss removal dramatically improved HMM quality — transitions
p50 dropped 60% and masked CNV rates fell from 40–71% to 12–33%. Evaluable FNR reached zero
for CRT/MDR1/PM2_PM3, confirming the sin loss was the root cause of noisy reconstructions.
Masked CRR is >> 1.55 for all genes (p10=1.58–1.62), confirming signal is present and the
CRR rescue threshold is not the bottleneck.

**Where predictions diverged:** GCH1 evaluable FNR=0.14 remains (244 FNs, all frac=0, CRR
p50=1.43 < 1.55 rescue threshold). The total FNR is still 0.12–0.33 because the CRR rescue
in `16_genome_cnv_caller` runs *after* the sanity check — sanity-failed samples get -1
before rescue fires, even when their CRR >> 1.55. Masked CNV pool persists for this reason,
not because of signal weakness.

**Next:** Experiment 43 — `18_genome_cnv_caller`: apply CRR rescue before the sanity check
so sanity-failed samples with CRR ≥ 1.55 are rescued rather than masked.

## Proposal history

**Original proposal (2026-05-13):** Smooth HMM emission input with a 5-bin centred rolling mean
via `04_gaussian_hmm`. Parameters: `hmm_smooth_window=5`, `hmm_self_transition=0.99` (unchanged),
reuse exp 37 checkpoint (no retrain). Rationale: self_transition raises had hit diminishing
returns; smoothing the emissions was the next lever to reduce HMM oscillation.

**Feedback round 1:** Remove the sin loss component from the error term — it is making
reconstructions really weird.

**What changed (round 1):** The proposal pivoted from an HMM-only fix to a VAE retrain without
sin loss. Rather than smoothing the noisy signal downstream, we fix the noise at its source.
The HMM smoothing (`04_gaussian_hmm`, `smooth_window=5`) was initially kept as a defensive
measure.

**Feedback round 2:** No rolling-mean HMM. Establish a clean baseline for 400bp.

**What changed (round 2):** Reverted HMM module to `03_gaussian_hmm` with no `smooth_window`.
The core change (sin loss removed, full VAE retrain) is unchanged. This is the minimal clean
baseline: only sin loss is removed; nothing else is added.

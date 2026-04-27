# Experiment Timeline

One-to-two line summary per experiment. Marks dead-ends (✗) and breakthroughs (★).

---

**Exp 01 — Baseline VAE + HMM**
Established 1000 bp bin size, latent dim 10, initial HMM. Severe overfitting; model learned CNV profiles as "normal." Poor GCH1 recall.

**Exp 02 — HMM stickiness tuning (self_transition 0.99→0.95)**
No retraining; loosened HMM to improve recall. Results not documented — led to exp 03.

**Exp 03 — Aggressive HMM tuning (self_transition 0.95→0.90)**
No retraining; pushed stickiness further. Results not documented — led to root-cause diagnosis in exp 04.

**Exp 04 — HMM state boundary fix + CRR fallback ★**
Identified CN=1/CN=2 boundary was at CRR=1.5, stranding FN signal. Fixed state mean to 1.5; added CRR fallback for short genes. PM2_PM3 recall recovered, but CRR fallback threshold 1.30 admitted excessive CRT false positives.

**Exp 05 — Bin-count-gated CRR fallback ✗**
Applied CRR fallback only to genes with <3 bins (physically unable to segment). Eliminated CRT/MDR1 FPs, but PM2_PM3 FNR shot to 0.98 — HMM still insufficient for this gene.

**Exp 06 — CNV sample downsampling (ratio 0.25) ★**
Full retrain. Weighted sampler reduced CNV prevalence during training, fixing biased VAE. Recall surged across all genes (PM2_PM3 FNR 0.98→0.20), but self_transition=0.75 — calibrated for weak VAE — became too loose (CRT 142 FPs).

**Exp 07 — Raise HMM self_transition 0.75→0.80**
No retrain. Modest FP reduction (142→118 CRT FPs) but borderline FP/TP CRR gap only 0.13; lever nearly exhausted.

**Exp 08 — Diagnostic baseline with segment callability metrics**
Eval-only upgrade (03_pf9_evaluation). Established callability=0.995, transitions p50=123. HMM oscillates rather than sustaining segments — characteristic of weak VAE signal.

**Exp 09 — Encoder dropout p=0.20 ★/✗**
Full retrain. GCH1 FNR improved 0.17→0.07, but PM2_PM3 collapsed 0.26→0.94. Dropout VAE produces cleaner but shorter reconstructions; HMM too sticky for short PM2_PM3 runs.

**Exp 10 — Lower HMM self_transition 0.80→0.75 ✗**
No retrain. PM2_PM3 barely changed (0.94→0.93). HMM tuning alone insufficient; FN CRR p50=1.43 confirms signal present but discarded.

**Exp 11 — Lower encoder dropout p=0.20→0.10 ✗**
Full retrain. All genes regressed simultaneously. p=0.10 shifts reconstruction error distribution; confounded comparison without re-tuning HMM.

**Exp 12 — Higher encoder dropout p=0.30 ★**
Full retrain. PM2_PM3 FNR dropped 0.93→0.53. Higher regularisation forces tighter normal-mode reconstruction, amplifying anomaly signal at inference. HMM remains the bottleneck.

**Exp 13 — Lower HMM self_transition 0.80→0.65**
No retrain. PM2_PM3 FNR dropped further 0.53→0.34. GCH1 FNR unchanged (FNs due to VAE signal, not HMM inertia). Separated two failure modes: PM2_PM3 is HMM-limited, GCH1 is VAE-limited.

**Exp 14 — Skip connections in VAE encoder (05_conv_vae) ★/✗**
Full retrain. PM2_PM3 FNR fell 0.34→0.08, but PPV collapsed everywhere (CRT 0.86→0.11). Skip connections amplified anomaly signal; self_transition=0.65 now far too loose.

**Exp 15 — Skip connections + self_transition 0.65→0.80**
No retrain. PM2_PM3 FNR rose to 0.24 (expected), but PPV did not recover (CRT=0.13). FPs are CNV-caller-level: many noise calls pass CRR gate 1.30, not a pure HMM problem.

**Exp 16 — Raise cnv_crr_amp_threshold 1.30→1.40 ✗**
No retrain. PPV unchanged for CRT/MDR1/PM2_PM3. Discovered the threshold only applies via the short-gene fallback path (<3 bins); long genes rely solely on HMM — threshold change was a no-op for them.

**Exp 17 — CRR gate on long genes (04_gene_cnv_caller) ★**
New caller adds CRR post-gate for all genes: HMM CN≥2 calls rejected if CRR<1.40. CRT PPV 0.13→0.91, MDR1 0.76→0.97. PM2_PM3 FNR=0.41 — gate cuts too many borderline TPs (CRR p50=1.28 among FNs).

**Exp 18 — Sin loss fine-tuning to sharpen integer signal ★/✗**
Fine-tuned from exp 14 checkpoint. GCH1 FNR improved 0.18→0.07; CRT/MDR1 near-perfect. PM2_PM3 collapsed 0.41→0.74 — sin loss sharpened segments into fragments; no single segment spans the gene, so strict spanning check never fires.

**Exp 19 — Fractional gene coverage caller + gate 1.40→1.35 (05_gene_cnv_caller) ★**
New caller replaces strict spanning check with ≥50% bin coverage rule. PM2_PM3 FNR 0.74→0.01. Both failure modes fixed: fragmented segments now aggregate, gate at 1.35 recovers borderline TPs. All four reference genes near-perfect.

**Exp 20 — Curriculum sampling (cnv ratio ramp 0.05→0.40) ✗**
Full retrain. GCH1 FNR improved only marginally (0.07→0.05); final ratio 0.40 over-exposed VAE to CNVs, inflating CRR for normals — PPV collapsed (CRT 0.91→0.71, 759 FPs).

**Exp 21 — Curriculum sampling, lower final ratio (0.05→0.30)**
Full retrain. GCH1 PPV recovered (0.71→0.91); FNR back to exp 19 baseline (0.07). 130 FNs have CRR p50=1.25 — signal genuinely weak, just below gate. Curriculum approach exhausted.

**Exp 22 — Synthetic normal augmentation (Poisson resampling) ✗**
Full retrain. GCH1 FNR 0.07→0.06 (trivial). FN CRR p50 actually fell 1.25→1.23. Augmented normals added noise rather than tightening the normal prior.

**Exp 23 — Genome-wide CNV calling + threshold 1.30→1.25**
New caller applied to all 5,318 genes. GCH1 FNR→0.04. CRT PPV fell (0.91→0.25); threshold 1.25 admits too many borderline calls. Diagnostic run establishing the genome-wide pipeline.

**Exp 24 — Threshold 1.25→1.20 + fix PM2_PM3 GFF lookup + logistic sharpness metric ✗**
No retrain. GCH1 FNR→0.02. CRT/MDR1 PPV collapsed further; HMM transitions p90=1167 (extreme jumpiness) newly observed. Threshold lowering backfired.

**Exp 25 — Raise HMM self_transition 0.80→0.90 to reduce jumpiness ✗**
No retrain. Transitions p90 fell only 1167→1024 (12%). PPV unchanged — FPs are systematic CRR elevation near threshold, not HMM micro-call noise.

**Exp 26 — Band segment filter: demote CRR [1.20, 1.25) calls lacking HMM segment support ★**
New caller: borderline CRR calls demoted to CN=1 unless ≥50% of gene bins covered by CN≥2 HMM segments. PPV recovered to exp 23 levels across all genes. Remaining FPs have CRR≥1.26 — artefact band extends above 1.25.

**Exp 27 — Widen artefact band filter 1.25→1.35 ★/✗**
Extended band upper limit. PPV gains exceeded predictions (GCH1 0.84→0.96, MDR1 0.71→0.91). But GCH1 FNR rose 0.04→0.11 — 223 genuine GCH1 amplifications at CRR p50=1.28 lack enough HMM segment coverage to pass the 0.50 threshold.

**Exp 28 — Two-tier band filter: core [1.20, 1.25)@0.50, extended [1.25, 1.35)@0.20 (proposed)**
Split-band caller relaxes segment coverage threshold in the extended zone (0.50→0.20) while keeping the strict gate in the core zone. Predicted to recover GCH1 FNR 0.11→0.05–0.07 with minor PPV cost.

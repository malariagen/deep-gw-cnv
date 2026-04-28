# models/

## Big picture — how the VAE signal works

The VAE produces a reconstruction of what "normal" read-count coverage should look like
for each sample, based on the overall coverage profile. At inference, CRR = input /
reconstruction: CNVs stick out because the reconstruction doesn't model them.

Two failure modes to always keep in mind:
- **Over-training**: the VAE memorises CNV profiles as normal → CNVs stop sticking out
  → HMM misses amplifications (high FNR).
- **Under-training**: reconstructions are too different from input everywhere → CRR is
  noisy → HMM produces many borderline calls (high FP).

CRR is a rough signal (gene/flank ratio over a 100 kb window). Do not over-optimise
HMM parameters to CRR percentiles — that is optimising a noisy proxy. When FP and FN
CRR distributions overlap, the problem is upstream (VAE reconstruction quality), not
a calling-layer problem. Address reconstruction quality first; tune the calling layer
after the signal is clean.

### Experiment design rules

**FN CRR p50 decision tree** — check this before proposing any experiment targeting FNs:
- p50 < 1.10 → genuinely weak signal; investigate data quality, not calling parameters.
- 1.10 ≤ p50 ≤ 1.40 → signal exists but soft; propose a **signal-boost training** experiment (see below). Do not lower the gate threshold — that trades FPs for FNs without fixing the root cause.
- p50 > 1.40 AND PPV has headroom → gate/threshold adjustment is valid.

**Signal-boosting strategies** (in priority order when FN CRR p50 is in 1.10–1.40):
1. Synthetic normal augmentation — Poisson-resample non-CNV samples during training so the VAE learns a tighter normal prior; true CNVs deviate more at inference.
2. Masked reconstruction loss — zero the loss on CNV gene regions so the VAE is not penalised for failing to reconstruct amplifications.
3. Architecture changes — more capacity or dropout for regularisation.

**HMM self-transition** is a cheap lever (no retraining). When PM2_PM3 FP=0 (precision headroom exists), try aggressive drops (e.g. 0.80→0.65) before retraining. Small steps (0.80→0.75) have historically been ineffective.

**`cnv_downsample_ratio`** — do not increase this to improve recall. Exposing the model to more CNV-positive samples shifts the learned normal baseline; the VAE starts reconstructing amplifications accurately, which kills the anomaly signal the HMM relies on. Only touch this parameter if you have a specific mechanistic reason it will help without shifting the baseline.

Architecture improvement directions (in rough priority order):
1. **Dropout** in encoder — prevents per-sample CNV memorisation (simplest, tried in exp 09)
2. **Masked training** — exclude CNV gene loci from reconstruction loss during training
3. **U-Net skip connections** in decoder — finer reconstruction in non-CNV regions
4. **Self-attention block** — captures long-range genomic dependencies

## Entry point

```bash
.venv/bin/python train.py path/to/experiment/config.yaml
```

Loads the config, dynamically imports the versioned components named in it, trains the model, then runs inference → HMM segmentation → CNV calling → evaluation in sequence.

## Modules

| Path | What it does |
|---|---|
| `train.py` | Reads `config.yaml`, loads versioned components, runs full pipeline |
| `architectures/` | `ConvVAE` — 1-D convolutional encoder/decoder (v01–v06) |
| `hmm/` | `run_hmm_all_samples`, `fit_hmm_sample` — Gaussian HMM segmentation (v01–v03) |
| `cnv/` | `run_cnv_calls`, `call_all_genes` — gene-level CNV calls from segments (v01–v11) |
| `evaluation/` | `run_evaluation` — metrics vs Pf9 GATK ground truth (v01–v03) |
| `training/dataset.py` | `ReadCountDataset` — reads NPY store, log2-normalises |
| `training/trainer.py` | `train_vae` — training loop with KL annealing + early stopping |
| `training/wrap_up.py` | Post-training pipeline: inference → HMM → CNV calling → evaluation |
| `experiments/` | One folder per experiment, each with a `config.yaml` and `run.sh` |

## Versioning

Each component family has its own folder with numbered implementation files (`01_*.py`, `02_*.py`, …). The `config.yaml` for each experiment declares which version to use:

```yaml
architecture: "01_conv_vae"
hmm:          "01_gaussian_hmm"
cnv:          "01_gene_cnv_caller"
evaluation:   "01_pf9_evaluation"
```

`train.py` loads these via `importlib` at runtime, so a new experiment can mix and match any existing versions just by editing the config. The `__init__.py` in each folder re-exports the latest version as a convenience for interactive use.

**Rule:** reuse an existing versioned file and change only the config parameters when only the parameters differ. Create a new numbered file only when the algorithm or interface itself changes.

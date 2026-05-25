# CLAUDE.md — orientation for every new session

Read this file completely before doing anything else. It explains what this project is,
how the experiment loop works, and the rules you must follow.

---

## What this project is

An autonomous ML research pipeline that discovers copy-number variants (CNVs) in
*P. falciparum* malaria whole-genome sequencing data. A VAE learns what "normal"
coverage looks like; deviations from the reconstruction signal amplifications/deletions.
The full pipeline: train VAE → run inference → HMM segmentation → CNV calling →
evaluate against Pf9 ground truth.

The pipeline runs unattended on a Mac mini. Claude proposes experiments by email, the
user replies AUTHORISE or sends feedback, and a launchd daemon runs the approved
experiment automatically.

---

## The autonomous experiment loop — READ THIS FIRST

This is the heart of the project. Every experiment goes through this loop.
**Never bypass it. Never run experiments directly.**

```
/propose-experiment  →  sends email proposal  →  user replies AUTHORISE
       ↑                                                   ↓
/propose-experiment  ←  experiment completes  ←  daemon runs run.sh
(proposes next exp)       (sends results email)   (tools/check_and_run.sh)
```

### Daemon state files (tools/)

| File | Purpose |
|---|---|
| `tools/.last_proposal_msgid` | Message-ID of the outstanding proposal email. Daemon checks for replies against this. **If this file exists, a proposal is live.** |
| `tools/.last_proposal_experiment` | Experiment number the daemon will run on AUTHORISE. |
| `tools/.proposal_thread_msgid` | Root Message-ID for the current feedback chain. Keeps all revisions in one email thread. |
| `tools/.proposal_body.txt` | Body of the most recent proposal (working file, not state). |
| `tools/daemon.log` | Daemon stdout/stderr. Check this when an experiment fails. |

The daemon (`tools/check_and_run.sh`) is installed as a launchd agent via
`bash tools/install_daemon.sh` and polls every 60 seconds.

### Claude commands (always use these, never ad-hoc scripts)

- **`/propose-experiment`** — analyse `evaluation.txt` from the last completed
  experiment, design the next experiment, create its folder, send the proposal email,
  and arm the daemon. Full instructions: `.claude/commands/propose-experiment.md`.
- **`/check-reply`** — check whether the user has replied AUTHORISE or sent feedback,
  and act on it. Full instructions: `.claude/commands/check-reply.md`.

### What happens on AUTHORISE

The daemon:
1. Sends a "starting" email.
2. Runs `bash models/experiments/N/run.sh >> tools/daemon.log 2>&1`.
3. On success: invokes `claude --print /propose-experiment` to propose the next experiment.
4. On failure: sends a failure email with a pointer to `daemon.log` and **stops the loop**.

On failure, Claude Code is opened, the bug is diagnosed, and once fixed,
`/propose-experiment` is run to re-arm the daemon. The loop restarts from there.

---

## How to handle a failed experiment

When an experiment fails (daemon sends a failure email):

1. **Read `tools/daemon.log`** to find the traceback.
2. **Fix the bug** — edit the relevant source file(s).
3. **If training already completed** (i.e. `checkpoint.pth` exists in the result dir):
   - Update `run.sh` to call `wrap_up.py` directly instead of `train.py`, so the
     daemon will skip training and resume from the checkpoint:
     ```bash
     .venv/bin/python -m training.wrap_up models/experiments/N/config.yaml
     ```
   - Delete any corrupt intermediate artifacts (e.g. a malformed `reconstructions.npy`).
   - Do NOT delete `checkpoint.pth`.
4. **If training had not yet started** (no `checkpoint.pth`): no changes to `run.sh`
   needed — it will retrain from scratch.
5. **Re-arm the daemon** by running `/propose-experiment`. It will detect that
   experiment N has no `evaluation.txt` (not yet complete), treat N−1 as the last
   completed, and propose N again (or a revised N).

**Never run `run.sh` or `train.py` manually. Never use Bash background tasks to run
experiments. Never send ad-hoc emails outside the loop.** Doing so breaks the
post-run proposal flow and leaves the user with no notification.

---

## Key files and directories

```
models/
  train.py                  Entry point. Loads config, runs full pipeline.
                            Skips training if checkpoint.pth exists (new behaviour);
                            skips inference too if reconstructions.npy exists.
  training/
    dataset.py              ReadCountDataset — reads NPY store, log2-normalises,
                            pads to ceil(n_bins/32)*32 for conv/FC layers.
    trainer.py              train_vae — training loop with KL annealing + early stopping.
    wrap_up.py              run_inference + post-training pipeline entry point.
                            Use as standalone to resume from checkpoint:
                            python -m training.wrap_up experiments/N/config.yaml
  architectures/            Versioned VAE definitions (01_conv_vae.py … 07_fc_vae.py)
  hmm/                      Versioned HMM segmenters  (01_ … 04_gaussian_hmm.py)
  cnv/                      Versioned CNV callers     (01_ … 18_genome_cnv_caller.py)
  evaluation/               Versioned evaluators      (01_ … 05_pf9_evaluation.py)
  experiments/
    N/
      config.yaml           Single source of truth: data paths, component versions,
                            all hyperparameters.
      run.sh                Runs the experiment. Must go through the daemon loop.
      README.md             Hypothesis, expected outcome, actual outcome (filled in
                            by /propose-experiment after completion).

tools/
  check_and_run.sh          Daemon script — polls email, runs experiments.
  check_reply.py            Checks Gmail for replies to a specific Message-ID.
  send_email.py             Sends email and optionally saves the Message-ID.
  install_daemon.sh         Installs the launchd agent (run once).
  daemon.log                All daemon and experiment stdout/stderr.

data/
  inputs/                   Read-count NPY stores (one per bin size / genome region).
  results/                  Per-experiment outputs, one folder per experiment.
    N_<name>/
      checkpoint.pth        Best model weights. Authoritative "training done" signal.
      reconstructions.npy   (n_samples, n_bins) denormalised reconstructions.
      segments.parquet       HMM segments.
      gene_calls.tsv        Gene-level CNV calls.
      evaluation.txt        Metrics vs ground truth. Existence = experiment complete.

assets/
  PlasmoDB-54_Pfalciparum3D7.gff   Gene annotations.
  20260313_full_cnv_data_pf9.tsv   Pf9 ground-truth CNV calls.
  Pf_9_samples_20260227.txt        Sample metadata.

diagnostics/app.py          Streamlit app for interactive inspection of any experiment.
.claude/commands/
  propose-experiment.md     Full instructions for /propose-experiment.
  check-reply.md            Full instructions for /check-reply.
```

---

## Coding rules

- Always run Python via `.venv/bin/python` — never system `python3` or `python`.
- Write simple code. Comment non-obvious logic. Update comments when editing code.
- Explicitly free large objects (`del df; gc.collect()`) once no longer needed.
- Never edit versioned files (numeric prefix `01_`, `02_`, …) in-place. Create a new
  numbered file for any algorithmic change. Parameter-only changes go in config.yaml.
- When a `run.sh` reuses a large artifact from a parent experiment unchanged
  (e.g. `segments.parquet`, `reconstructions.npy`), symlink it (`ln -sf`), don't copy.
- The CNV caller must use a single global CRR rescue threshold — never per-gene
  thresholds. The caller is genome-wide; gene-specific thresholds do not generalise.
- Never run experiments directly (via Bash or background tasks). All experiment
  execution goes through the daemon email loop. See "How to handle a failed experiment"
  above for the correct recovery procedure.

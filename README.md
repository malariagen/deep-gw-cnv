# AI autoresearch for copy number variation (CNV) discovery in malaria field samples using variational autoencoders

In March, Karpathy (legendary AI researcher) released [autoresearch](https://github.com/karpathy/autoresearch), where an AI agent continuously proposes and performs 5-minute ML experiments to improve a mini LLM. The idea was that you can go to sleep and by the time you wake up, your ML model has already explored 100 different architectures/hyperparameters/datasets without you having to constantly decide for yourself. Does this "autoresearch" strategy translate to complex biological problems? 

Malaria kills 600,000 people a year. Amplifications of certain genes in malaria have already made some drugs be dropped from the first-line treatment regimen in multiple countries. I'll talk through how I'm debating with AIs over email and getting my Mac mini to work for me as close to 24/7 as possible, to discover novel CNVs by "autoresearching" a variational autoencoder-based pipeline on WGS data. 

![screenshot](assets/screenshots/Screenshot%202026-04-08%20at%2017.59.23.png)
![screenshot](assets/screenshots/Screenshot%202026-04-08%20at%2018.01.55.png)

## Layout

```
data/
  inputs/       read-count NPY stores
  results/      per-experiment outputs (one folder per experiment)
  setup/        scripts to extract read counts from BAMs/CRAMs

assets/         sample manifests, BAM/CRAM path lists, reference files

models/
  train.py          entry point — runs a full experiment from a config
  architectures/    versioned VAE definitions  (01_conv_vae.py, …)
  hmm/              versioned HMM segmenters   (01_gaussian_hmm.py, …)
  cnv/              versioned CNV callers       (01_gene_cnv_caller.py, …)
  evaluation/       versioned evaluators        (01_pf9_evaluation.py, …)
  training/         dataset loader, trainer, inference (non-versioned)
  experiments/      one self-contained folder per experiment

diagnostics/    Streamlit app for interactive sample inspection
```

## How it works

```mermaid
flowchart LR
    subgraph Experiment["models/experiments/N/"]
        direction TB
        cfg["config.yaml<br/>arch · hmm · cnv · eval · paths"]
        vc["Versioned components"]
        cfg -->|"selects versions & params"| vc
    end

    subgraph Pipeline["Execution"]
        direction TB
        runsh["run.sh"]
        out["data/results/N/"]
        runsh -->|"inference → HMM → CNV → eval"| out
    end

    subgraph Loop["Autonomous proposal loop"]
        direction TB
        claude["Claude Code<br/>/propose-experiment"]
        mail["📧 Proposal email"]
        you(["You"])
        daemon["Daemon · launchd · 60 s"]
        claude --> mail --> you
        you -->|"AUTHORISE"| daemon
        you -->|"feedback"| claude
    end

    diag["Streamlit app<br/>diagnostics/app.py"]

    vc --> runsh
    out -->|"evaluation.txt"| claude
    claude -->|"creates N+1 config"| cfg
    daemon -->|"runs experiment"| runsh
    cfg -. "resolves paths & versions" .-> diag
    out -.-> diag
```

## Experiment proposal workflow

Claude analyses the latest `evaluation.txt`, proposes the next experiment, creates the folder, and emails a summary. Reply "AUTHORISE" to run it on the Mac mini; reply with feedback to get a revised proposal.

**First-time setup** (install the background polling daemon):
```bash
bash tools/install_daemon.sh
```

**To propose the next experiment** (invoke from Claude Code):
```
/propose-experiment
```
Claude sets up the experiment folder, writes a README, and emails a ≤100-line summary.

**The daemon** (`tools/check_and_run.sh`, running via launchd every 60 s) checks for a reply:
- `AUTHORISE` → runs the experiment automatically
- Anything else → flags that feedback is waiting; open Claude Code and run `/check-reply`

**Privacy:** `check_reply.py` searches only by the exact Message-ID of the proposal email. It never lists or reads any other email. Reply body content is never written to disk or logs.

**VSCode vs Desktop:** Either works. The daemon runs independently of which editor is open.

## Diagnostics

```bash
cd diagnostics
streamlit run app.py
```

Select an experiment from the dropdown — the app loads that experiment's config to resolve data paths, component versions, and all calling parameters automatically.

## Setup

See [data/setup/](data/setup/) for scripts that extract read counts from BAMs/CRAMs and convert to NPY.

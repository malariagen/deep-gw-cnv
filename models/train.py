"""
Central training entry point.

Usage:
    python train.py path/to/config.yaml

All paths in the config are resolved relative to the config file's directory,
so experiments stay self-contained regardless of where you invoke this script.
"""

import argparse
import inspect
import os
import sys

import importlib

import numpy as np
import pandas as pd
import yaml
import torch
from torch.utils.data import DataLoader, WeightedRandomSampler

# Allow running this file from any working directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from training.dataset import AugmentedNormalDataset, ReadCountDataset
from training.trainer import train_vae
from training.wrap_up import run_inference


def get_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def _make_downsampled_loader(ds, cfg, resolve, downsample_ratio, batch_size, num_workers):
    """Return a DataLoader that down-weights CNV-positive samples.

    CNV-positive samples are drawn at `downsample_ratio` relative to normal
    samples (e.g. 0.25 → sampled 4× less often). This prevents the VAE from
    learning amplified profiles as the default state in datasets where CNVs
    are common (e.g. PM2 in Pf9).
    """
    gt_path = cfg.get("pf9_gt_path")
    if not gt_path:
        raise ValueError("cnv_downsample_ratio requires pf9_gt_path to be set in config")

    gt = pd.read_csv(resolve(gt_path), sep="\t", index_col="Sample")

    # Any sample with at least one positive final amplification or deletion call
    call_cols = [c for c in gt.columns if c.endswith("_final_amplification_call")
                 or c.endswith("_final_deletion_call")]
    cnv_positive = set(gt.index[(gt[call_cols] == 1).any(axis=1)])

    sample_ids = ds.sample_ids
    weights = np.where(
        np.isin(sample_ids, list(cnv_positive)),
        downsample_ratio,
        1.0,
    ).astype(np.float64)

    n_cnv  = int((weights < 1.0).sum())
    n_norm = int((weights == 1.0).sum())
    print(f"CNV downsampling: {n_cnv} CNV-positive samples (ratio={downsample_ratio}), "
          f"{n_norm} normal samples", flush=True)

    sampler = WeightedRandomSampler(weights, num_samples=len(weights), replacement=True)
    return DataLoader(ds, batch_size=batch_size, sampler=sampler, num_workers=num_workers)


def _make_curriculum_loader_fn(ds, cfg, resolve, batch_size, num_workers):
    """Return a per-epoch DataLoader factory for curriculum-weighted sampling.

    Linearly ramps cnv_downsample_ratio from cnv_downsample_ratio_initial to
    cnv_downsample_ratio_final over cnv_downsample_warmup_epochs epochs.
    CNV-positive samples start underrepresented so the VAE learns normal profiles
    first, then receive increasing weight so the VAE learns CNV profiles too.
    """
    ratio_initial = cfg["cnv_downsample_ratio_initial"]
    ratio_final   = cfg["cnv_downsample_ratio_final"]
    warmup        = cfg["cnv_downsample_warmup_epochs"]

    gt_path = cfg.get("pf9_gt_path")
    if not gt_path:
        raise ValueError("curriculum sampling requires pf9_gt_path to be set in config")

    gt = pd.read_csv(resolve(gt_path), sep="\t", index_col="Sample")
    call_cols = [c for c in gt.columns if c.endswith("_final_amplification_call")
                 or c.endswith("_final_deletion_call")]
    cnv_positive = set(gt.index[(gt[call_cols] == 1).any(axis=1)])
    is_cnv = np.isin(ds.sample_ids, list(cnv_positive))

    n_cnv  = int(is_cnv.sum())
    n_norm = int((~is_cnv).sum())
    print(f"Curriculum sampler: {n_cnv} CNV-positive, {n_norm} normal | "
          f"ratio {ratio_initial:.2f}→{ratio_final:.2f} over {warmup} epochs", flush=True)

    def make_loader(epoch):
        t     = min(1.0, epoch / max(warmup, 1))
        ratio = ratio_initial + t * (ratio_final - ratio_initial)
        weights = np.where(is_cnv, ratio, 1.0).astype(np.float64)
        if epoch % 10 == 0:
            print(f"  [curriculum] epoch {epoch}: cnv_ratio={ratio:.3f}", flush=True)
        sampler = WeightedRandomSampler(weights, num_samples=len(weights), replacement=True)
        return DataLoader(ds, batch_size=batch_size, sampler=sampler, num_workers=num_workers)

    return make_loader


def main():
    parser = argparse.ArgumentParser(description="Train a model from a config file.")
    parser.add_argument("config", help="Path to experiment config.yaml")
    args = parser.parse_args()

    config_dir = os.path.dirname(os.path.abspath(args.config))
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    def resolve(path):
        return path if os.path.isabs(path) else os.path.join(config_dir, path)

    store_path = resolve(cfg["store_path"])
    out_dir    = resolve(cfg["out_dir"])

    # Load versioned components named in config
    ConvVAE           = importlib.import_module(f"architectures.{cfg['architecture']}").ConvVAE
    run_hmm_all_samples = importlib.import_module(f"hmm.{cfg['hmm']}").run_hmm_all_samples
    run_cnv_calls       = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
    run_evaluation      = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

    lsf_jobid = os.getenv("LSB_JOBID")
    if lsf_jobid:
        print(f"LSF job: {lsf_jobid}\n{'='*50}", flush=True)

    device = get_device()
    print(f"Device: {device}", flush=True)

    ds = ReadCountDataset(store_path, normalise=cfg["normalise"])

    # Optionally wrap with synthetic normal augmentation.
    # Poisson-resamples each normal sample per __getitem__ call, so every epoch
    # sees a fresh draw. CNV-positive samples are passed through unchanged.
    if cfg.get("aug_normal_poisson"):
        gt_path = cfg.get("pf9_gt_path")
        if not gt_path:
            raise ValueError("aug_normal_poisson requires pf9_gt_path in config")
        gt = pd.read_csv(resolve(gt_path), sep="\t", index_col="Sample")
        call_cols  = [c for c in gt.columns if c.endswith("_final_amplification_call")
                      or c.endswith("_final_deletion_call")]
        cnv_ids    = set(gt.index[(gt[call_cols] == 1).any(axis=1)])
        is_cnv     = np.isin(ds.sample_ids, list(cnv_ids))
        scale_range = cfg.get("aug_normal_depth_scale")  # e.g. [0.85, 1.15] or None
        n_aug  = int((~is_cnv).sum())
        n_cnv  = int(is_cnv.sum())
        print(f"Normal augmentation: {n_aug} normal samples will be Poisson-resampled "
              f"({n_cnv} CNV-positive samples unchanged). "
              f"depth_scale={scale_range}", flush=True)
        ds = AugmentedNormalDataset(ds, is_cnv, depth_scale_range=scale_range)

    num_workers = 8 if device.type == "cuda" else 0

    make_loader_fn   = None
    downsample_ratio = cfg.get("cnv_downsample_ratio")
    if cfg.get("cnv_downsample_ratio_initial") is not None:
        # Curriculum: ramp CNV weight from initial to final over warmup epochs.
        make_loader_fn = _make_curriculum_loader_fn(ds, cfg, resolve, cfg["batch_size"], num_workers)
        dl = make_loader_fn(0)
    elif downsample_ratio is not None:
        dl = _make_downsampled_loader(ds, cfg, resolve, downsample_ratio, cfg["batch_size"], num_workers)
    else:
        dl = DataLoader(ds, batch_size=cfg["batch_size"], shuffle=True, num_workers=num_workers)

    # Architectures >= 06 accept n_bins; older ones take only latent_dim.
    if 'n_bins' in inspect.signature(ConvVAE.__init__).parameters:
        model = ConvVAE(n_bins=ds.n_bins_raw, latent_dim=cfg["latent_dim"]).to(device)
    else:
        model = ConvVAE(latent_dim=cfg["latent_dim"]).to(device)

    # Optionally load a pre-trained checkpoint before fine-tuning
    pretrained = cfg.get("pretrained_checkpoint")
    if pretrained:
        pretrained = pretrained if os.path.isabs(pretrained) else os.path.join(config_dir, pretrained)
        model.load_state_dict(torch.load(pretrained, map_location=device, weights_only=True))
        print(f"Loaded pretrained checkpoint: {pretrained}", flush=True)

    optimiser = torch.optim.Adam(
        model.parameters(),
        lr=cfg["lr"],
        weight_decay=cfg["weight_decay"],
    )

    print(f"samples={len(ds)} | latent_dim={cfg['latent_dim']}", flush=True)

    os.makedirs(out_dir, exist_ok=True)
    checkpoint_path = os.path.join(out_dir, "checkpoint.pth")
    log_path        = os.path.join(out_dir, "training_log.json")

    train_vae(
        model, dl, optimiser,
        epochs                 = cfg["epochs"],
        max_beta               = cfg["max_beta"],
        warmup_epochs          = cfg["warmup_epochs"],
        patience               = cfg["patience"],
        device                 = device,
        model_save_path        = checkpoint_path,
        log_path               = log_path,
        sin_loss_max_weight    = cfg.get("sin_loss_max_weight", 0.0),
        sin_loss_warmup_epochs = cfg.get("sin_loss_warmup_epochs", 0),
        make_loader_fn         = make_loader_fn,
    )

    if os.path.exists(checkpoint_path):
        model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
        print("Loaded best checkpoint for inference.", flush=True)

    run_inference(model, ds, device, out_dir, batch_size=cfg["batch_size"])

    print("Fitting HMM segments...", flush=True)
    run_hmm_all_samples(store_path, out_dir, cfg)

    print("Calling gene CNVs...", flush=True)
    run_cnv_calls(store_path, out_dir, cfg)

    if cfg.get("pf9_gt_path"):
        cfg_resolved = dict(cfg)
        cfg_resolved["pf9_gt_path"] = resolve(cfg["pf9_gt_path"])
        if cfg.get("pf9_meta_path"):
            cfg_resolved["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])
        run_evaluation(out_dir, cfg_resolved)
    else:
        print("Skipping evaluation (pf9_gt_path not set in config).", flush=True)


if __name__ == "__main__":
    main()

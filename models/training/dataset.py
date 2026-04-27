import math
import numpy as np
import torch
import torch.nn.functional as F
import os
from torch.utils.data import Dataset


class ReadCountDataset(Dataset):
    def __init__(self, npy_dir: str, normalise: bool = True):
        # npy_dir should contain counts.npy and sample_ids.npy.
        # mmap_mode='r' avoids loading the full array into RAM at startup — important
        # for high-resolution datasets (e.g. 400 bp core: ~11 GB).
        self.counts     = np.load(os.path.join(npy_dir, "counts.npy"), mmap_mode='r')
        self.sample_ids = np.load(os.path.join(npy_dir, "sample_ids.npy"), allow_pickle=True)
        self.normalise  = normalise
        self.n_bins_raw    = self.counts.shape[1]
        # Pad to the next multiple of 32 so 5× stride-2 conv blocks divide evenly.
        self.n_bins_padded = math.ceil(self.n_bins_raw / 32) * 32

    def __len__(self):
        return len(self.sample_ids)

    def __getitem__(self, idx):
        counts = self.counts[idx, :].astype(np.float32)
        if self.normalise:
            counts = np.log2(counts + 1)
        t = torch.from_numpy(counts)
        t = F.pad(t, (0, self.n_bins_padded - self.n_bins_raw))
        return t


class AugmentedNormalDataset(Dataset):
    """Wraps ReadCountDataset and applies Poisson resampling to normal samples.

    Normal samples are resampled from Poisson(expected_counts) each epoch,
    producing diverse coverage profiles that tighten the VAE's normal prior.
    CNV-positive samples are passed through unchanged.

    The goal is higher CRR at inference for true CNVs: a VAE trained on richer
    normal variation reconstructs unseen CNV profiles less faithfully, amplifying
    the anomaly signal without changing the detection gate.
    """

    def __init__(self, base_dataset: ReadCountDataset, is_cnv: np.ndarray,
                 depth_scale_range=None):
        """
        base_dataset       : ReadCountDataset (normalise=True produces log2 inputs).
        is_cnv             : bool array of shape (N,) — True for CNV-positive samples.
        depth_scale_range  : optional (lo, hi) tuple; depth is scaled by U[lo, hi]
                             before Poisson sampling to simulate library-size variation.
        """
        self.ds                = base_dataset
        self.is_cnv            = is_cnv
        self.depth_scale_range = depth_scale_range
        self.sample_ids        = base_dataset.sample_ids

    def __len__(self):
        return len(self.ds)

    def __getitem__(self, idx):
        item = self.ds[idx]  # shape (n_bins_padded,), log2(count+1) normalised

        if self.is_cnv[idx]:
            return item

        # Convert log2 space back to expected raw counts, then Poisson-resample.
        # Padded bins are 0 → 2^0 - 1 = 0 → Poisson(0) = 0 always; safe.
        raw = (2.0 ** item - 1.0).clamp(min=0.0)

        if self.depth_scale_range is not None:
            lo, hi = self.depth_scale_range
            scale  = torch.empty(1).uniform_(lo, hi).item()
            raw    = raw * scale

        # torch.poisson requires float input
        augmented = torch.poisson(raw)
        return torch.log2(augmented + 1.0)

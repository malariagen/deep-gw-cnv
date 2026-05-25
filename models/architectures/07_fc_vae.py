"""
07_fc_vae — Fully-connected (MLP) VAE
======================================

Replaces the 1-D convolutional blocks of 06_conv_vae with dense layers,
removing the spatial-locality inductive bias entirely.

With convolutional layers the encoder can learn "elevated counts around
a specific region is normal in African samples" by specialising local
filters on those bins.  A fully-connected encoder applies the same
global weight matrix to the entire coverage vector, making it structurally
harder to memorise population-specific regional patterns like the prevalent
African GCH1 amplification.  The hope is that the bottleneck forces higher
reconstruction residuals for any profile that deviates from the global mean,
including GCH1-amplified African samples that currently have frac=0.

Architecture:
  Encoder:  Linear(n_bins_padded, 2048) → BN → ReLU → Dropout(0.30)
            Linear(2048, 512)           → BN → ReLU → Dropout(0.30)
            Linear(512, 128)            → BN → ReLU
            mu     head: Linear(128, latent_dim)
            logvar head: Linear(128, latent_dim)

  Decoder:  Linear(latent_dim, 128)     → BN → ReLU
            Linear(128, 512)            → BN → ReLU
            Linear(512, 2048)           → BN → ReLU
            Linear(2048, n_bins_padded) (no activation — raw reconstruction)

n_bins_padded = ceil(n_bins / 32) * 32, matching the padding applied by
ReadCountDataset so the first linear layer dimension is consistent.

ConvVAE is exported as an alias so train.py's generic import works unchanged.
"""

import math
import torch
import torch.nn as nn


class FCEncoder(nn.Module):
    def __init__(self, n_bins_padded: int, latent_dim: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_bins_padded, 2048),
            nn.BatchNorm1d(2048),
            nn.ReLU(),
            nn.Dropout(p=0.30),
            nn.Linear(2048, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(p=0.30),
            nn.Linear(512, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
        )
        self.mu     = nn.Linear(128, latent_dim)
        self.logvar = nn.Linear(128, latent_dim)

    def forward(self, x: torch.Tensor):
        h      = self.net(x)
        mu     = self.mu(h)
        logvar = torch.clamp(self.logvar(h), -10, 10)
        return mu, logvar


class FCDecoder(nn.Module):
    def __init__(self, n_bins_padded: int, latent_dim: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Linear(128, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Linear(512, 2048),
            nn.BatchNorm1d(2048),
            nn.ReLU(),
            nn.Linear(2048, n_bins_padded),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.net(z)


class FCVAE(nn.Module):
    def __init__(self, n_bins: int, latent_dim: int):
        super().__init__()
        # Match ReadCountDataset's padding so linear layer dimensions are consistent.
        n_bins_padded = math.ceil(n_bins / 32) * 32
        self.n_bins = n_bins
        self.enc = FCEncoder(n_bins_padded, latent_dim)
        self.dec = FCDecoder(n_bins_padded, latent_dim)

    def forward(self, x: torch.Tensor) -> dict:
        mu, logvar = self.enc(x)
        z          = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        # Trim padding so reconstruction shape matches raw counts (same as conv VAE).
        return {"recon": self.dec(z)[:, :self.n_bins], "z": (mu, logvar)}


# Alias for train.py compatibility
ConvVAE = FCVAE

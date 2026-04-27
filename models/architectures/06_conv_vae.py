"""
06_conv_vae — 1-D Convolutional VAE with dynamic n_bins
=========================================================

Identical to 05_conv_vae (residual skip connections in encoder) except that
the number of genome bins is a constructor argument rather than a module-level
constant.  This lets the same architecture train on any bin resolution (e.g.
400 bp or 1000 bp) without code changes: training code derives n_bins from
the dataset shape and passes it in.

  ConvVAE(n_bins=51986, latent_dim=10)   # 400 bp core genome
  ConvVAE(n_bins=20814, latent_dim=10)   # 1000 bp core genome (05_conv_vae equivalent)

n_bins_padded = ceil(n_bins / 32) * 32 is computed internally so the 5×
stride-2 encoder always produces an integer spatial dimension after downsampling.

Architecture (unchanged from 05_conv_vae):
  Encoder:  5× ResConvBlock, channels 1→32→64→128→256→256,
            main path: Conv1d(k=7, stride=2) → BN → ReLU → Dropout(p=0.30),
            shortcut:  Conv1d(k=1, stride=2) → BN.
  Decoder:  5× ConvTranspose1d blocks, no dropout.
"""

import math
import torch
import torch.nn as nn


# ---------------------------------------------------------------------------
# Encoder block with residual skip connection
# ---------------------------------------------------------------------------

class ResConvBlock(nn.Module):
    """Stride-2 conv block with a 1×1 strided shortcut to preserve spatial detail."""

    def __init__(self, in_ch: int, out_ch: int, stride: int = 2):
        super().__init__()
        self.main = nn.Sequential(
            nn.Conv1d(in_ch, out_ch, kernel_size=7, stride=stride, padding=3),
            nn.BatchNorm1d(out_ch),
            nn.ReLU(),
            nn.Dropout(p=0.30),
        )
        # 1×1 conv to match output spatial size and channel count for the residual add
        self.shortcut = nn.Sequential(
            nn.Conv1d(in_ch, out_ch, kernel_size=1, stride=stride),
            nn.BatchNorm1d(out_ch),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.main(x) + self.shortcut(x)


# ---------------------------------------------------------------------------
# Encoder / Decoder
# ---------------------------------------------------------------------------

class ConvEncoder(nn.Module):
    def __init__(self, n_bins_padded: int, latent_dim: int):
        super().__init__()
        self.blocks = nn.Sequential(
            ResConvBlock(1,   32,  stride=2),
            ResConvBlock(32,  64,  stride=2),
            ResConvBlock(64,  128, stride=2),
            ResConvBlock(128, 256, stride=2),
            ResConvBlock(256, 256, stride=2),
        )
        self.flat_dim = 256 * (n_bins_padded // 32)
        self.mu     = nn.Linear(self.flat_dim, latent_dim)
        self.logvar = nn.Linear(self.flat_dim, latent_dim)

    def forward(self, x: torch.Tensor):
        h      = self.blocks(x.unsqueeze(1)).flatten(1)
        mu     = self.mu(h)
        logvar = torch.clamp(self.logvar(h), -10, 10)
        return mu, logvar


class ConvDecoder(nn.Module):
    def __init__(self, n_bins: int, n_bins_padded: int, latent_dim: int):
        super().__init__()
        self.n_bins        = n_bins
        self.n_bins_padded = n_bins_padded
        self.flat_dim      = 256 * (n_bins_padded // 32)
        self.proj          = nn.Linear(latent_dim, self.flat_dim)
        self.deconv        = nn.Sequential(
            self._block(256, 256, stride=2),
            self._block(256, 128, stride=2),
            self._block(128,  64, stride=2),
            self._block( 64,  32, stride=2),
            self._block( 32,   1, stride=2),
        )

    @staticmethod
    def _block(in_ch: int, out_ch: int, stride: int) -> nn.Sequential:
        return nn.Sequential(
            nn.ConvTranspose1d(in_ch, out_ch, kernel_size=7, stride=stride,
                               padding=3, output_padding=1),
            nn.BatchNorm1d(out_ch),
            nn.ReLU(),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        h     = self.proj(z).view(z.size(0), 256, self.n_bins_padded // 32)
        recon = self.deconv(h).squeeze(1)
        return recon[:, :self.n_bins]


# ---------------------------------------------------------------------------
# VAE
# ---------------------------------------------------------------------------

class ConvVAE(nn.Module):
    def __init__(self, n_bins: int, latent_dim: int):
        super().__init__()
        n_bins_padded = math.ceil(n_bins / 32) * 32
        self.enc = ConvEncoder(n_bins_padded, latent_dim)
        self.dec = ConvDecoder(n_bins, n_bins_padded, latent_dim)

    def forward(self, x: torch.Tensor) -> dict:
        mu, logvar = self.enc(x)
        z          = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return {"recon": self.dec(z), "z": (mu, logvar)}

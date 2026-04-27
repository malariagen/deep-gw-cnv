#!/usr/bin/env python3
"""
Rebin Pf9-53973-samples-100bp-npy to 400bp and filter to core genome.

Run from the project root:
  .venv/bin/python data/setup/100bp_to_400bp_core.py

Outputs:
  data/inputs/Pf9-53973-samples-400bp-npy       (all 400bp bins)
  data/inputs/Pf9-53973-samples-400bp-core-npy  (400bp bins overlapping core regions)
"""

import gc
import time
import numpy as np
from pathlib import Path

SRC_DIR  = Path("data/inputs/Pf9-53973-samples-100bp-npy")
DST_400  = Path("data/inputs/Pf9-53973-samples-400bp-npy")
DST_CORE = Path("data/inputs/Pf9-53973-samples-400bp-core-npy")
CORE_BED = Path("assets/core-genome-onebased.bed")

# Samples processed per iteration. Each chunk uses ~(CHUNK * n_bins_100 * 8) bytes
# of RAM for the uint64 intermediate (~370 MB at CHUNK=200).
CHUNK = 200


def build_400bp_bins(contigs_100):
    """
    Group every 4 consecutive 100bp bins within each chromosome into one 400bp bin.

    Returns:
        contigs_400  - structured array (same dtype as contigs_100)
        group_starts - 1D int array of the start index in contigs_100 for each 400bp
                       group, suitable for np.add.reduceat
    """
    new_contigs  = []
    group_starts = []

    i = 0
    n = len(contigs_100)
    while i < n:
        chrom = contigs_100[i]["chrom"]
        # Advance j to the end of this chromosome's bins.
        j = i
        while j < n and contigs_100[j]["chrom"] == chrom:
            j += 1
        # Slice bins i..j-1 into groups of 4 (last group may be partial).
        k = i
        while k < j:
            g_end = min(k + 4, j)
            group_starts.append(k)
            start_bp = int(contigs_100[k]["start"])
            end_bp   = int(contigs_100[g_end - 1]["end"])
            new_contigs.append((chrom, start_bp, end_bp))
            k += 4
        i = j

    contigs_400  = np.array(new_contigs, dtype=contigs_100.dtype)
    group_starts = np.array(group_starts, dtype=np.intp)
    return contigs_400, group_starts


def core_overlap_mask(contigs_400, bed_path):
    """
    Return a boolean mask (len == len(contigs_400)) that is True for every
    400bp bin that shares at least 1bp overlap with a core-genome region.

    Both BED and contigs use 1-based inclusive coordinates.
    """
    core = {}
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            core.setdefault(chrom, []).append((start, end))
    for ch in core:
        core[ch].sort()

    n    = len(contigs_400)
    mask = np.zeros(n, dtype=bool)
    for i in range(n):
        chrom     = contigs_400[i]["chrom"]
        intervals = core.get(chrom)
        if not intervals:
            continue
        bin_s = int(contigs_400[i]["start"])
        bin_e = int(contigs_400[i]["end"])
        for ivl_s, ivl_e in intervals:
            if ivl_s > bin_e:
                break  # sorted; no further interval can overlap
            if ivl_e >= bin_s:
                mask[i] = True
                break
    return mask


def main():
    print("Loading 100bp contigs and sample IDs…")
    contigs_100 = np.load(SRC_DIR / "contigs.npy", allow_pickle=True)
    sample_ids  = np.load(SRC_DIR / "sample_ids.npy", allow_pickle=True)
    n_samples   = len(sample_ids)
    n_bins_100  = len(contigs_100)
    print(f"  {n_samples} samples, {n_bins_100} 100bp bins")

    print("Building 400bp bin layout…")
    contigs_400, group_starts = build_400bp_bins(contigs_100)
    n_bins_400 = len(contigs_400)
    print(f"  {n_bins_400} 400bp bins")

    print(f"Computing core genome overlap mask from {CORE_BED} (1-based)…")
    mask    = core_overlap_mask(contigs_400, CORE_BED)
    n_core  = int(mask.sum())
    core_idx = np.where(mask)[0]
    contigs_core = contigs_400[mask]
    print(f"  {n_core} 400bp bins overlap core regions")

    DST_400.mkdir(parents=True, exist_ok=True)
    DST_CORE.mkdir(parents=True, exist_ok=True)

    np.save(DST_400  / "contigs.npy",    contigs_400)
    np.save(DST_400  / "sample_ids.npy", sample_ids)
    np.save(DST_CORE / "contigs.npy",    contigs_core)
    np.save(DST_CORE / "sample_ids.npy", sample_ids)
    print("Saved contigs.npy and sample_ids.npy for both outputs.")

    print("Opening source counts (mmap, read-only)…")
    counts_100 = np.load(str(SRC_DIR / "counts.npy"), mmap_mode="r")

    print("Creating output count memmaps…")
    out_400  = np.lib.format.open_memmap(
        str(DST_400 / "counts.npy"),  mode="w+", dtype="uint32",
        shape=(n_samples, n_bins_400),
    )
    out_core = np.lib.format.open_memmap(
        str(DST_CORE / "counts.npy"), mode="w+", dtype="uint32",
        shape=(n_samples, n_core),
    )

    print(f"Processing {n_samples} samples in chunks of {CHUNK}…")
    t0 = time.time()
    for chunk_start in range(0, n_samples, CHUNK):
        chunk_end  = min(chunk_start + CHUNK, n_samples)

        # Cast to uint64 before summing to avoid overflow of 4× uint32 values.
        chunk = counts_100[chunk_start:chunk_end].astype(np.uint64)

        # Rebin: sum consecutive 100bp groups using reduceat along the bin axis.
        # Each group_starts[g] is the index of the first 100bp bin in 400bp group g;
        # reduceat sums bins [group_starts[g], group_starts[g+1]) automatically.
        chunk_400 = np.add.reduceat(chunk, group_starts, axis=1).astype(np.uint32)
        del chunk
        gc.collect()

        out_400[chunk_start:chunk_end]  = chunk_400
        out_core[chunk_start:chunk_end] = chunk_400[:, core_idx]
        del chunk_400
        gc.collect()

        elapsed = time.time() - t0
        rate    = chunk_end / max(elapsed, 1e-3)
        eta     = (n_samples - chunk_end) / rate
        print(f"  {chunk_end}/{n_samples} ({chunk_end / n_samples * 100:.1f}%)"
              f" — {elapsed:.0f}s elapsed, ETA {eta:.0f}s")

    # Flush memmaps to disk.
    del out_400, out_core
    gc.collect()

    print("Done.")
    print(f"  400bp full : {DST_400}")
    print(f"  400bp core : {DST_CORE}")


if __name__ == "__main__":
    main()

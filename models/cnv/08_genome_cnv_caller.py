"""
Genome-wide gene CNV caller — version 08.

Extends version 07 by fixing a normalization artefact in the logistic sharpness metric.
In v07, both x (rank) and y (copy ratio) were normalised to [0, 1] before fitting the
sigmoid.  This allowed a flat-but-noisy sample — where a few bins are slightly elevated
above an otherwise uniform background — to produce a high slope: after y-normalisation,
the "peak" bins reach y=1, mimicking a clean bimodal boundary.

The fix: anchor the sigmoid center at the expected gene/flank boundary position rather
than letting it float.  If there are N_flank flanking bins and N_gene gene bins in the
local window, the transition in the *sorted* distribution should occur at rank
x ≈ N_flank / (N_flank + N_gene): flank bins are expected to occupy the lower-ranked
portion, gene bins the upper.  By fixing the center, only a transition that occurs
*at the right place* in the sorted sequence earns a high slope.

The fixed center is computed from the bin counts *before* NaN removal so that it
reflects the structural layout of the locus, not the noise pattern of any one sample.

`logistic_slope` is written to both output files:
  genome_wide_amplifications.tsv — one value per amplified (sample, gene) call.
  gene_calls.tsv                 — one value per sample per reference gene
                                   (all samples, including CN=1, for diagnostics).

All other calling logic is identical to version 07.
"""

import gc
import os
import re
import time
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


# Only nuclear chromosomes; skip apicoplast (API) and mitochondrion (MIT).
NUCLEAR_CHROMS = {f"Pf3D7_{i:02d}_v3" for i in range(1, 15)}

# Reference genes kept for backward-compatible gene_calls.tsv / evaluation.
# PM2_PM3 spans two separate GFF genes (PMII / PF3D7_1408000 and PMIII / PF3D7_1408100);
# both are tagged call_id "PM2_PM3" and per-sample CN/CRR is the max of the two.
REFERENCE_GENES = [
    {"call_id": "MDR1",    "contig": "Pf3D7_05_v3", "start": 955955,  "end": 963095},
    {"call_id": "CRT",     "contig": "Pf3D7_07_v3", "start": 402385,  "end": 406341},
    {"call_id": "GCH1",    "contig": "Pf3D7_12_v3", "start": 974226,  "end": 976097},
    {"call_id": "PM2_PM3", "contig": "Pf3D7_14_v3", "start": 292244,  "end": 295261},  # PMII
    {"call_id": "PM2_PM3", "contig": "Pf3D7_14_v3", "start": 296683,  "end": 299101},  # PMIII
]

# Build (contig, start, end) → call_id for O(1) lookup while iterating genes.
_REF_LOOKUP = {(g["contig"], g["start"], g["end"]): g["call_id"] for g in REFERENCE_GENES}


def _parse_gff_genes(gff_path):
    """Return list of dicts {gene_id, gene_name, contig, start, end} for nuclear coding genes."""
    attr_re = re.compile(r'(?:^|;)(\w+)=([^;]+)')
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'protein_coding_gene':
                continue
            contig = parts[0]
            if contig not in NUCLEAR_CHROMS:
                continue
            attrs     = dict(attr_re.findall(parts[8]))
            gene_id   = attrs.get('ID', '')
            gene_name = attrs.get('Name', gene_id)
            genes.append({
                'gene_id':   gene_id,
                'gene_name': gene_name,
                'contig':    contig,
                'start':     int(parts[3]),
                'end':       int(parts[4]),
            })
    return genes


def _precompute_sanity_ok(segments, sample_ids, sid_to_idx, chroms, starts,
                           min_confidence, min_cn1_proportion):
    """Return dict {chrom: bool_array[n_samples]} — True where chromosome sanity passes.

    A sample passes the sanity check for a chromosome if >= min_cn1_proportion of
    that chromosome's bins are covered by high-confidence CN=1 HMM segments.
    Uses binary search (searchsorted) to avoid O(n_bins * n_segs) inner loops.
    """
    n = len(sample_ids)
    sanity_ok = {}

    for chrom in NUCLEAR_CHROMS:
        chrom_mask = chroms == chrom
        if not chrom_mask.any():
            sanity_ok[chrom] = np.zeros(n, dtype=bool)
            continue

        s_chrom = np.sort(starts[chrom_mask])
        n_chrom = len(s_chrom)

        hc1 = segments[
            (segments["chrom"] == chrom) &
            (segments["cn"] == 1) &
            (segments["confidence"] >= min_confidence)
        ]

        ok = np.zeros(n, dtype=bool)
        for sid, grp in hc1.groupby("sample_id"):
            idx = sid_to_idx.get(sid)
            if idx is None:
                continue
            covered = np.zeros(n_chrom, dtype=bool)
            for x0, x1 in zip(grp["x0"].values, grp["x1"].values):
                lo = int(np.searchsorted(s_chrom, x0, side='left'))
                hi = int(np.searchsorted(s_chrom, x1, side='left'))
                covered[lo:hi] = True
            if covered.sum() / n_chrom >= min_cn1_proportion:
                ok[idx] = True
        sanity_ok[chrom] = ok

    return sanity_ok


def _logistic_sharpness(crr_bins, n_flank_expected, n_total_expected):
    """Fit a logistic sigmoid to sorted bin copy ratios with a fixed center.

    The center is anchored at n_flank_expected / n_total_expected — the expected
    fraction of the sorted distribution occupied by flanking (CN=1) bins.  This
    prevents the normalization artefact in v07: a flat sample with a few noisy peaks
    can no longer fool the metric by placing the inflection near the tail of the
    sorted y-distribution, because the center is fixed at the structural boundary,
    not at wherever the best-fit inflection happens to fall.

    Returns nan when the fit fails or fewer than 6 finite bins are available.
    """
    n = len(crr_bins)
    if n < 6:
        return np.nan

    sorted_v = np.sort(crr_bins)
    v_min, v_max = sorted_v[0], sorted_v[-1]
    if v_max - v_min < 1e-6:
        # No variation — unimodal; sharpness is undefined / not applicable.
        return 0.0

    x = np.linspace(0.0, 1.0, n)
    y = (sorted_v - v_min) / (v_max - v_min)

    # Fix center at the expected flank/gene boundary; only slope is free.
    center = n_flank_expected / n_total_expected

    def _sigmoid(x, slope):
        return 1.0 / (1.0 + np.exp(-np.clip(slope * (x - center), -50.0, 50.0)))

    try:
        popt, _ = curve_fit(
            _sigmoid, x, y,
            p0=[5.0],
            bounds=([0.0], [500.0]),
            maxfev=500,
        )
        return float(popt[0])
    except Exception:
        return np.nan


def run_cnv_calls(store_path, out_dir, cfg):
    """Compute genome-wide gene CNV calls and write outputs."""
    min_cn1_proportion = cfg["cnv_min_cn1_proportion"]
    min_confidence     = cfg["cnv_min_confidence"]
    flank_padding      = cfg["cnv_flank_padding"]
    crr_amp_threshold  = cfg["cnv_crr_amp_threshold"]
    gff_path           = cfg["gff_path"]

    contigs    = pd.DataFrame(np.load(os.path.join(store_path, "contigs.npy"), allow_pickle=True))
    counts     = np.load(os.path.join(store_path, "counts.npy"))
    recons     = np.load(os.path.join(out_dir, "reconstructions.npy"))
    sample_ids = np.load(os.path.join(out_dir, "sample_ids.npy"), allow_pickle=True)
    segments   = pd.read_parquet(os.path.join(out_dir, "segments.parquet"))

    chroms = contigs["chrom"].values
    starts = contigs["start"].values.astype(float)
    del contigs
    gc.collect()

    copy_ratios = counts.astype(float) / (recons.astype(float) + 1e-6)
    del counts, recons
    gc.collect()

    n = len(sample_ids)
    sid_to_idx = {sid: i for i, sid in enumerate(sample_ids)}

    print("Precomputing per-chromosome sanity checks...", flush=True)
    sanity_ok = _precompute_sanity_ok(
        segments, sample_ids, sid_to_idx,
        chroms, starts, min_confidence, min_cn1_proportion,
    )
    del segments
    gc.collect()

    print(f"Parsing genes from {gff_path}...", flush=True)
    genes = _parse_gff_genes(gff_path)
    print(f"  {len(genes):,} protein-coding genes on nuclear chromosomes", flush=True)

    genes_by_chrom = defaultdict(list)
    for gene in genes:
        genes_by_chrom[gene["contig"]].append(gene)

    print(f"Calling genome-wide CNVs for {n:,} samples...", flush=True)
    t0 = time.time()

    amp_rows      = []   # → genome_wide_amplifications.tsv
    ref_gene_rows = []   # → gene_calls.tsv (4 reference genes only)

    n_genes_done = 0
    for chrom in sorted(genes_by_chrom.keys()):
        chrom_genes = genes_by_chrom[chrom]
        chrom_mask  = chroms == chrom
        s_chrom     = starts[chrom_mask]
        cr_chrom    = copy_ratios[:, chrom_mask]   # [n_samples, n_chrom_bins]
        ok          = sanity_ok.get(chrom, np.zeros(n, dtype=bool))

        for gene in chrom_genes:
            g_start   = gene["start"]
            g_end     = gene["end"]
            gene_id   = gene["gene_id"]
            gene_name = gene["gene_name"]

            gene_mask_l   = (s_chrom >= g_start) & (s_chrom <= g_end)
            flank_mask_l  = (s_chrom < g_start - flank_padding) | (s_chrom > g_end + flank_padding)
            # Window for sharpness metric: gene + padding on both sides.
            window_mask_l = (s_chrom >= g_start - flank_padding) & (s_chrom <= g_end + flank_padding)

            if not gene_mask_l.any() or not flank_mask_l.any():
                continue

            # Bin counts for fixed-center logistic (structural layout, not per-sample).
            # gene_mask_l is a subset of window_mask_l; flank bins = window - gene.
            n_gene_in_window  = int(gene_mask_l.sum())
            n_total_in_window = int(window_mask_l.sum())
            n_flank_in_window = n_total_in_window - n_gene_in_window

            mean_gene  = np.nanmean(cr_chrom[:, gene_mask_l],  axis=1)
            mean_flank = np.nanmean(cr_chrom[:, flank_mask_l], axis=1)
            crr_all    = np.where(mean_flank > 0, mean_gene / mean_flank, np.nan)

            # Vectorised call: missing → -1, CRR >= threshold → 2, else → 1.
            finite_mask = np.isfinite(crr_all)
            amp_mask    = ok & finite_mask & (crr_all >= crr_amp_threshold)

            amp_indices = np.where(amp_mask)[0]
            for idx in amp_indices:
                crr_window_bins = cr_chrom[idx, window_mask_l]
                finite_bins     = crr_window_bins[np.isfinite(crr_window_bins)]
                slope           = _logistic_sharpness(
                    finite_bins, n_flank_in_window, n_total_in_window,
                )
                amp_rows.append({
                    "sample_id":      sample_ids[idx],
                    "gene_id":        gene_id,
                    "gene_name":      gene_name,
                    "chrom":          chrom,
                    "start":          g_start,
                    "end":            g_end,
                    "crr":            round(float(crr_all[idx]), 4),
                    "logistic_slope": round(float(slope), 4) if np.isfinite(slope) else None,
                })

            call_id = _REF_LOOKUP.get((chrom, g_start, g_end))
            if call_id is not None:
                # Pre-compute sharpness for all samples on this reference gene.
                slopes = np.full(n, np.nan)
                for idx in range(n):
                    crr_window_bins = cr_chrom[idx, window_mask_l]
                    finite_bins     = crr_window_bins[np.isfinite(crr_window_bins)]
                    slopes[idx]     = _logistic_sharpness(
                        finite_bins, n_flank_in_window, n_total_in_window,
                    )

                for idx in range(n):
                    crr_val = float(crr_all[idx]) if finite_mask[idx] else None
                    if not ok[idx]:
                        cn = -1
                    elif crr_val is not None and crr_val >= crr_amp_threshold:
                        cn = 2
                    else:
                        cn = 1
                    slope_val = float(slopes[idx]) if np.isfinite(slopes[idx]) else None
                    ref_gene_rows.append({
                        "sample_id":      sample_ids[idx],
                        "call_id":        call_id,
                        "cn":             cn,
                        "crr":            crr_val,
                        "logistic_slope": round(slope_val, 4) if slope_val is not None else None,
                    })

            n_genes_done += 1
            if n_genes_done % 500 == 0:
                print(f"  {n_genes_done}/{len(genes)} genes | elapsed {time.time() - t0:.0f}s",
                      flush=True)

        del cr_chrom
        gc.collect()

    print(f"Done | {len(genes):,} genes processed | elapsed {time.time() - t0:.0f}s", flush=True)

    del copy_ratios, sanity_ok
    gc.collect()

    # Write genome-wide amplifications master file (long format, CN>=2 only).
    if amp_rows:
        amp_df = pd.DataFrame(amp_rows)
        amp_df = amp_df.sort_values(["chrom", "start", "sample_id"]).reset_index(drop=True)
    else:
        amp_df = pd.DataFrame(
            columns=["sample_id", "gene_id", "gene_name", "chrom", "start", "end",
                     "crr", "logistic_slope"]
        )
    amp_path = os.path.join(out_dir, "genome_wide_amplifications.tsv")
    amp_df.to_csv(amp_path, sep="\t", index=False)
    n_unique_samples = amp_df["sample_id"].nunique() if len(amp_df) else 0
    n_unique_genes   = amp_df["gene_id"].nunique()   if len(amp_df) else 0
    print(
        f"Saved genome-wide amplifications: {len(amp_df):,} calls | "
        f"{n_unique_samples:,} samples | {n_unique_genes:,} genes → {amp_path}",
        flush=True,
    )
    del amp_rows, amp_df
    gc.collect()

    # Write gene_calls.tsv (wide format for the 4 reference genes, for evaluation).
    ref_df = pd.DataFrame(ref_gene_rows)
    # Aggregate duplicate call_ids (e.g. PM2_PM3 spans two GFF genes: PMII + PMIII).
    # CN: max (CN=2 if either sub-gene is amplified). CRR and slope: max of the two.
    ref_agg = (
        ref_df.groupby(["sample_id", "call_id"])
              .agg(
                  cn=("cn", "max"),
                  crr=("crr", lambda x: x.max(skipna=True)),
                  logistic_slope=("logistic_slope", lambda x: x.max(skipna=True)),
              )
              .reset_index()
    )
    cn_wide    = ref_agg.pivot(index="sample_id", columns="call_id", values="cn")
    crr_wide   = ref_agg.pivot(index="sample_id", columns="call_id", values="crr")
    slope_wide = ref_agg.pivot(index="sample_id", columns="call_id", values="logistic_slope")
    cn_wide.columns.name    = None
    crr_wide.columns.name   = None
    slope_wide.columns.name = None

    # Ensure every reference gene has a column even if coordinates didn't match the GFF.
    ref_ids = list(dict.fromkeys(g["call_id"] for g in REFERENCE_GENES))  # unique, ordered
    missing = [g for g in ref_ids if g not in cn_wide.columns]
    if missing:
        print(f"WARNING: reference gene(s) {missing} not matched by GFF coordinates — "
              "columns will be NaN in gene_calls.tsv", flush=True)
        for g in missing:
            cn_wide[g]    = pd.NA
            crr_wide[g]   = np.nan
            slope_wide[g] = np.nan

    crr_wide.columns   = [f"crr_{g}"            for g in crr_wide.columns]
    slope_wide.columns = [f"logistic_slope_{g}"  for g in slope_wide.columns]
    gene_calls_wide = pd.concat([cn_wide, crr_wide, slope_wide], axis=1).reset_index()
    gene_path = os.path.join(out_dir, "gene_calls.tsv")
    gene_calls_wide.to_csv(gene_path, sep="\t", index=False)
    print(f"Saved reference gene calls ({len(gene_calls_wide):,} samples) → {gene_path}", flush=True)
    del ref_gene_rows, ref_df, ref_agg, cn_wide, crr_wide, slope_wide, gene_calls_wide
    gc.collect()

"""
Genome-wide gene CNV caller — version 11.

Extends v10 with a two-tier band filter.  Exp 27 (v10, band_upper=1.35,
cn2_threshold=0.50) raised GCH1 FNR from 0.04 to 0.11 because genuine
amplifications at CRR 1.25–1.35 have only partial HMM CN>=2 coverage —
the signal is weaker at that amplitude and the HMM doesn't fully cover the
gene body.  A single global threshold of 0.50 is too strict for this zone.

v11 splits the filtered band into two zones with independent thresholds:

  Core zone    [crr_amp_threshold, cnv_crr_band_core_upper):
      requires CN>=2 fraction >= cnv_crr_band_cn2_threshold (e.g. 0.50).
      Artefact-band FPs in this range typically have fraction ≈ 0, so 0.50
      is a clean discriminator (proven in exp 26).

  Extended zone [cnv_crr_band_core_upper, cnv_crr_band_upper):
      requires CN>=2 fraction >= cnv_crr_band_ext_cn2_threshold (e.g. 0.20).
      Genuine amplifications in this range have weaker signal; a lower
      threshold recovers TPs with partial HMM coverage while still blocking
      pure artefacts (fraction ≈ 0).

Calls with CRR >= cnv_crr_band_upper are admitted unconditionally (unchanged).

New config keys
--------------------------
cnv_crr_band_core_upper        — CRR boundary between core and extended zones.
                                 Defaults to cnv_crr_band_upper (single-zone,
                                 backward-compatible with v10 behaviour).
cnv_crr_band_ext_cn2_threshold — CN>=2 fraction threshold for the extended zone.
                                 Defaults to cnv_crr_band_cn2_threshold.
"""

import gc
import os
import re
import time
from collections import defaultdict

import numpy as np
import pandas as pd


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


def _segment_cn2_fraction(cn2_intervals, gene_bin_starts):
    """Return fraction of gene bins covered by CN>=2 HMM segments.

    cn2_intervals: list of (x0, x1) tuples for CN>=2 segments of this (sample, chrom).
                   x0 is inclusive, x1 is exclusive (same convention as searchsorted usage).
    gene_bin_starts: sorted array of bin start positions within the gene.

    Returns 0.0 if no CN>=2 segments exist for this sample on this chrom.
    Returns nan if gene_bin_starts is empty.
    """
    n = len(gene_bin_starts)
    if n == 0:
        return np.nan
    if not cn2_intervals:
        return 0.0

    covered = np.zeros(n, dtype=bool)
    for x0, x1 in cn2_intervals:
        lo = int(np.searchsorted(gene_bin_starts, x0, side='left'))
        hi = int(np.searchsorted(gene_bin_starts, x1, side='left'))
        covered[lo:hi] = True

    return covered.sum() / n


def run_cnv_calls(store_path, out_dir, cfg):
    """Compute genome-wide gene CNV calls and write outputs."""
    min_cn1_proportion = cfg["cnv_min_cn1_proportion"]
    min_confidence     = cfg["cnv_min_confidence"]
    flank_padding      = cfg["cnv_flank_padding"]
    crr_amp_threshold  = cfg["cnv_crr_amp_threshold"]
    gff_path           = cfg["gff_path"]

    # Two-tier band segment filter (v11):
    #   Core zone    [crr_amp_threshold, band_core_upper): threshold = band_cn2_threshold
    #   Extended zone [band_core_upper, band_upper):       threshold = band_ext_cn2_threshold
    # Calls with CRR >= band_upper are admitted unconditionally.
    band_upper             = float(cfg.get("cnv_crr_band_upper", float("inf")))
    band_cn2_threshold     = float(cfg.get("cnv_crr_band_cn2_threshold", 0.0))
    # Defaults make v11 backward-compatible with v10 (single-zone behaviour).
    band_core_upper        = float(cfg.get("cnv_crr_band_core_upper", band_upper))
    band_ext_cn2_threshold = float(cfg.get("cnv_crr_band_ext_cn2_threshold", band_cn2_threshold))
    band_filter_on         = band_upper > crr_amp_threshold and band_cn2_threshold > 0.0

    contigs    = pd.DataFrame(np.load(os.path.join(store_path, "contigs.npy"), allow_pickle=True))
    counts     = np.load(os.path.join(store_path, "counts.npy"))
    recons     = np.load(os.path.join(out_dir, "reconstructions.npy"))
    sample_ids = np.load(os.path.join(out_dir, "sample_ids.npy"), allow_pickle=True)
    # Segments are kept in scope (not deleted after sanity check) so the band filter
    # can look up per-sample CN>=2 coverage at each gene.
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
    # Do NOT delete segments here — needed for band filter below.
    gc.collect()

    print(f"Parsing genes from {gff_path}...", flush=True)
    genes = _parse_gff_genes(gff_path)
    print(f"  {len(genes):,} protein-coding genes on nuclear chromosomes", flush=True)

    if band_filter_on:
        if band_core_upper < band_upper:
            print(
                f"Band segment filter ON (two-tier): "
                f"core [{crr_amp_threshold}, {band_core_upper}) threshold={band_cn2_threshold}, "
                f"extended [{band_core_upper}, {band_upper}) threshold={band_ext_cn2_threshold}",
                flush=True,
            )
        else:
            print(
                f"Band segment filter ON: calls in CRR [{crr_amp_threshold}, {band_upper}) "
                f"require CN>=2 gene coverage fraction >= {band_cn2_threshold}",
                flush=True,
            )

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

        # Build per-sample CN>=2 segment lookup for this chromosome.
        # Used by the band filter and by reference-gene diagnostics.
        chrom_cn2 = segments[(segments["chrom"] == chrom) & (segments["cn"] >= 2)]
        cn2_by_sid = {
            sid: list(zip(grp["x0"].values, grp["x1"].values))
            for sid, grp in chrom_cn2.groupby("sample_id")
        }
        del chrom_cn2
        gc.collect()

        for gene in chrom_genes:
            g_start   = gene["start"]
            g_end     = gene["end"]
            gene_id   = gene["gene_id"]
            gene_name = gene["gene_name"]

            gene_mask_l   = (s_chrom >= g_start) & (s_chrom <= g_end)
            flank_mask_l  = (s_chrom < g_start - flank_padding) | (s_chrom > g_end + flank_padding)

            if not gene_mask_l.any() or not flank_mask_l.any():
                continue

            gene_bin_starts = s_chrom[gene_mask_l]   # sorted; used for segment coverage

            mean_gene  = np.nanmean(cr_chrom[:, gene_mask_l],  axis=1)
            mean_flank = np.nanmean(cr_chrom[:, flank_mask_l], axis=1)
            crr_all    = np.where(mean_flank > 0, mean_gene / mean_flank, np.nan)

            finite_mask = np.isfinite(crr_all)
            amp_mask    = ok & finite_mask & (crr_all >= crr_amp_threshold)

            amp_indices = np.where(amp_mask)[0]
            for idx in amp_indices:
                sid      = sample_ids[idx]
                cn2_frac = _segment_cn2_fraction(cn2_by_sid.get(sid, []), gene_bin_starts)
                # Two-tier band filter: core zone uses band_cn2_threshold, extended zone
                # uses the more lenient band_ext_cn2_threshold.
                if band_filter_on and crr_all[idx] < band_upper:
                    thr = band_cn2_threshold if crr_all[idx] < band_core_upper else band_ext_cn2_threshold
                    if not np.isfinite(cn2_frac) or cn2_frac < thr:
                        continue
                amp_rows.append({
                    "sample_id":            sid,
                    "gene_id":              gene_id,
                    "gene_name":            gene_name,
                    "chrom":                chrom,
                    "start":                g_start,
                    "end":                  g_end,
                    "crr":                  round(float(crr_all[idx]), 4),
                    "segment_cn2_fraction": round(float(cn2_frac), 4) if np.isfinite(cn2_frac) else None,
                })

            call_id = _REF_LOOKUP.get((chrom, g_start, g_end))
            if call_id is not None:
                # Compute segment_cn2_fraction for all samples on this reference gene
                # (including CN=1 samples) so the column is available for diagnostics.
                cn2_fracs = np.full(n, np.nan)
                for idx in range(n):
                    sid = sample_ids[idx]
                    cn2_fracs[idx] = _segment_cn2_fraction(
                        cn2_by_sid.get(sid, []), gene_bin_starts,
                    )

                for idx in range(n):
                    crr_val  = float(crr_all[idx]) if finite_mask[idx] else None
                    frac_val = float(cn2_fracs[idx]) if np.isfinite(cn2_fracs[idx]) else None
                    if not ok[idx]:
                        cn = -1
                    elif crr_val is not None and crr_val >= crr_amp_threshold:
                        # Two-tier band filter for reference genes.
                        if band_filter_on and crr_val < band_upper:
                            thr = band_cn2_threshold if crr_val < band_core_upper else band_ext_cn2_threshold
                            cn = 2 if (frac_val is not None and frac_val >= thr) else 1
                        else:
                            cn = 2
                    else:
                        cn = 1
                    ref_gene_rows.append({
                        "sample_id":            sample_ids[idx],
                        "call_id":              call_id,
                        "cn":                   cn,
                        "crr":                  crr_val,
                        "segment_cn2_fraction": round(frac_val, 4) if frac_val is not None else None,
                    })

            n_genes_done += 1
            if n_genes_done % 500 == 0:
                print(f"  {n_genes_done}/{len(genes)} genes | elapsed {time.time() - t0:.0f}s",
                      flush=True)

        del cr_chrom, cn2_by_sid
        gc.collect()

    print(f"Done | {len(genes):,} genes processed | elapsed {time.time() - t0:.0f}s", flush=True)

    del copy_ratios, sanity_ok, segments
    gc.collect()

    # Write genome-wide amplifications master file (long format, CN>=2 only).
    if amp_rows:
        amp_df = pd.DataFrame(amp_rows)
        amp_df = amp_df.sort_values(["chrom", "start", "sample_id"]).reset_index(drop=True)
    else:
        amp_df = pd.DataFrame(
            columns=["sample_id", "gene_id", "gene_name", "chrom", "start", "end",
                     "crr", "segment_cn2_fraction"]
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
    # CN: max (CN=2 if either sub-gene is amplified). CRR and fraction: max of the two.
    ref_agg = (
        ref_df.groupby(["sample_id", "call_id"])
              .agg(
                  cn=("cn", "max"),
                  crr=("crr", lambda x: x.max(skipna=True)),
                  segment_cn2_fraction=("segment_cn2_fraction", lambda x: x.max(skipna=True)),
              )
              .reset_index()
    )
    cn_wide    = ref_agg.pivot(index="sample_id", columns="call_id", values="cn")
    crr_wide   = ref_agg.pivot(index="sample_id", columns="call_id", values="crr")
    frac_wide  = ref_agg.pivot(index="sample_id", columns="call_id", values="segment_cn2_fraction")
    cn_wide.columns.name   = None
    crr_wide.columns.name  = None
    frac_wide.columns.name = None

    # Ensure every reference gene has a column even if coordinates didn't match the GFF.
    ref_ids = list(dict.fromkeys(g["call_id"] for g in REFERENCE_GENES))  # unique, ordered
    missing = [g for g in ref_ids if g not in cn_wide.columns]
    if missing:
        print(f"WARNING: reference gene(s) {missing} not matched by GFF coordinates — "
              "columns will be NaN in gene_calls.tsv", flush=True)
        for g in missing:
            cn_wide[g]   = pd.NA
            crr_wide[g]  = np.nan
            frac_wide[g] = np.nan

    crr_wide.columns  = [f"crr_{g}"                   for g in crr_wide.columns]
    frac_wide.columns = [f"segment_cn2_fraction_{g}"  for g in frac_wide.columns]
    gene_calls_wide = pd.concat([cn_wide, crr_wide, frac_wide], axis=1).reset_index()
    gene_path = os.path.join(out_dir, "gene_calls.tsv")
    gene_calls_wide.to_csv(gene_path, sep="\t", index=False)
    print(f"Saved reference gene calls ({len(gene_calls_wide):,} samples) → {gene_path}", flush=True)
    del ref_gene_rows, ref_df, ref_agg, cn_wide, crr_wide, frac_wide, gene_calls_wide
    gc.collect()

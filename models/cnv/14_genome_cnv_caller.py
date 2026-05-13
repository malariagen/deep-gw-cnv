"""
Genome-wide gene CNV caller — version 14.

Extends v13 with two clarifications:

1. Confidence threshold lowered: cnv_min_confidence is now applied at 0.25 (down from 0.50
   in v13). More CN>=2 segments become eligible for coverage counting; this should recover
   FNs that had real signal but low-confidence HMM segments.

2. CN>=3 transitions within gene body handled explicitly: the elevated-CN interval set
   passed to _elevated_cn_fraction includes ALL segments with cn >= 2 (i.e. CN=2, CN=3,
   CN=4, ...). If the HMM transitions from CN=2 to CN=3 mid-gene, both segments contribute
   to the coverage fraction. The union of all elevated-CN intervals is computed so that
   multi-state amplifications are treated identically to uniform CN=2 amplifications.
   The output column is still named 'segment_cn2_fraction' (shorthand for 'fraction of
   gene bins in CN>=2 state') to preserve compatibility with the evaluation module.

In v14, the primary call criterion remains:
    >= cnv_min_gene_coverage_fraction of the gene's bins are covered by elevated-CN
    (cn >= 2, confidence >= cnv_min_confidence) HMM segments.

CRR is still computed and stored as a diagnostic column; it does not influence calls.
"""

import gc
import os
import re
import time
from collections import defaultdict

import numpy as np
import pandas as pd


NUCLEAR_CHROMS = {f"Pf3D7_{i:02d}_v3" for i in range(1, 15)}

REFERENCE_GENES = [
    {"call_id": "MDR1",    "contig": "Pf3D7_05_v3", "start": 955955,  "end": 963095},
    {"call_id": "CRT",     "contig": "Pf3D7_07_v3", "start": 402385,  "end": 406341},
    {"call_id": "GCH1",    "contig": "Pf3D7_12_v3", "start": 974226,  "end": 976097},
    {"call_id": "PM2_PM3", "contig": "Pf3D7_14_v3", "start": 292244,  "end": 295261},  # PMII
    {"call_id": "PM2_PM3", "contig": "Pf3D7_14_v3", "start": 296683,  "end": 299101},  # PMIII
]

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

    A sample passes if >= min_cn1_proportion of that chromosome's bins are covered by
    high-confidence CN=1 HMM segments.
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


def _elevated_cn_fraction(elevated_intervals, gene_bin_starts):
    """Return fraction of gene bins covered by elevated-CN (cn>=2) HMM segments.

    elevated_intervals: list of (x0, x1) tuples for all cn>=2 segments, regardless
        of whether they are CN=2, CN=3, CN=4, etc. The union of intervals is used, so
        a CN=2→CN=3 transition within the gene body is handled correctly: both contribute
        to the fraction.
    gene_bin_starts: sorted array of bin start positions within the gene.
    Returns nan if gene_bin_starts is empty, 0.0 if no elevated-CN intervals.
    """
    n = len(gene_bin_starts)
    if n == 0:
        return np.nan
    if not elevated_intervals:
        return 0.0

    covered = np.zeros(n, dtype=bool)
    for x0, x1 in elevated_intervals:
        lo = int(np.searchsorted(gene_bin_starts, x0, side='left'))
        hi = int(np.searchsorted(gene_bin_starts, x1, side='left'))
        covered[lo:hi] = True

    return covered.sum() / n


def run_cnv_calls(store_path, out_dir, cfg):
    """Compute genome-wide gene CNV calls and write outputs."""
    min_cn1_proportion       = cfg["cnv_min_cn1_proportion"]
    min_confidence           = cfg["cnv_min_confidence"]
    min_gene_coverage_frac   = cfg["cnv_min_gene_coverage_fraction"]
    flank_padding            = cfg["cnv_flank_padding"]
    gff_path                 = cfg["gff_path"]
    # Reuse HMM low-coverage threshold: bins where counts OR reconstruction fall below this
    # are excluded from the gene/flank copy-ratio means (same v12 fix; prevents CRR collapse).
    low_cov_threshold        = int(cfg.get("hmm_low_cov_threshold", 10))

    contigs    = pd.DataFrame(np.load(os.path.join(store_path, "contigs.npy"), allow_pickle=True))
    counts     = np.load(os.path.join(store_path, "counts.npy"),   mmap_mode='r')
    recons     = np.load(os.path.join(out_dir, "reconstructions.npy"), mmap_mode='r')
    sample_ids = np.load(os.path.join(out_dir, "sample_ids.npy"), allow_pickle=True)
    segments   = pd.read_parquet(os.path.join(out_dir, "segments.parquet"))

    chroms = contigs["chrom"].values
    starts = contigs["start"].values.astype(float)
    del contigs
    gc.collect()

    counts_f = counts.astype(np.float32);  del counts;  gc.collect()
    recons_f = recons.astype(np.float32);  del recons;  gc.collect()
    low_cov_mask = (counts_f < low_cov_threshold) | (recons_f < low_cov_threshold)
    recons_f += np.float32(1e-6)
    np.divide(counts_f, recons_f, out=counts_f)
    copy_ratios = counts_f
    del recons_f
    copy_ratios[low_cov_mask] = np.nan
    del low_cov_mask
    gc.collect()

    n = len(sample_ids)
    sid_to_idx = {sid: i for i, sid in enumerate(sample_ids)}

    print("Precomputing per-chromosome sanity checks...", flush=True)
    sanity_ok = _precompute_sanity_ok(
        segments, sample_ids, sid_to_idx,
        chroms, starts, min_confidence, min_cn1_proportion,
    )
    gc.collect()

    print(f"Parsing genes from {gff_path}...", flush=True)
    genes = _parse_gff_genes(gff_path)
    print(f"  {len(genes):,} protein-coding genes on nuclear chromosomes", flush=True)
    print(
        f"HMM-primary calling: CN>=2 gene-body coverage >= {min_gene_coverage_frac:.2f} "
        f"(confidence >= {min_confidence}); CRR stored as diagnostic only.",
        flush=True,
    )

    genes_by_chrom = defaultdict(list)
    for gene in genes:
        genes_by_chrom[gene["contig"]].append(gene)

    print(f"Calling genome-wide CNVs for {n:,} samples...", flush=True)
    t0 = time.time()

    amp_rows      = []
    ref_gene_rows = []

    n_genes_done = 0
    for chrom in sorted(genes_by_chrom.keys()):
        chrom_genes = genes_by_chrom[chrom]
        chrom_mask  = chroms == chrom
        s_chrom     = starts[chrom_mask]
        cr_chrom    = copy_ratios[:, chrom_mask]
        ok          = sanity_ok.get(chrom, np.zeros(n, dtype=bool))

        # All elevated-CN (cn>=2) segments for this chromosome, confidence-filtered.
        # Includes CN=3, CN=4, etc.: a CN=2→CN=3 transition within a gene is handled
        # by unioning all intervals in _elevated_cn_fraction.
        chrom_elevated = segments[
            (segments["chrom"] == chrom) &
            (segments["cn"] >= 2) &
            (segments["confidence"] >= min_confidence)
        ]
        cn_elevated_by_sid = {
            sid: list(zip(grp["x0"].values, grp["x1"].values))
            for sid, grp in chrom_elevated.groupby("sample_id")
        }
        del chrom_elevated
        gc.collect()

        for gene in chrom_genes:
            g_start   = gene["start"]
            g_end     = gene["end"]
            gene_id   = gene["gene_id"]
            gene_name = gene["gene_name"]

            gene_mask_l  = (s_chrom >= g_start) & (s_chrom <= g_end)
            flank_mask_l = (s_chrom < g_start - flank_padding) | (s_chrom > g_end + flank_padding)

            if not gene_mask_l.any() or not flank_mask_l.any():
                continue

            gene_bin_starts = s_chrom[gene_mask_l]

            # CRR: diagnostic only; not used for the call decision.
            mean_gene  = np.nanmean(cr_chrom[:, gene_mask_l],  axis=1)
            mean_flank = np.nanmean(cr_chrom[:, flank_mask_l], axis=1)
            crr_all    = np.where(mean_flank > 0, mean_gene / mean_flank, np.nan)
            finite_crr = np.isfinite(crr_all)

            # Genome-wide calls: only iterate samples with elevated-CN segments —
            # vastly more efficient than scanning all n samples.
            for sid, elevated_intervals in cn_elevated_by_sid.items():
                idx = sid_to_idx.get(sid)
                if idx is None or not ok[idx]:
                    continue
                elev_frac = _elevated_cn_fraction(elevated_intervals, gene_bin_starts)
                if not np.isfinite(elev_frac) or elev_frac < min_gene_coverage_frac:
                    continue
                crr_val = round(float(crr_all[idx]), 4) if finite_crr[idx] else None
                amp_rows.append({
                    "sample_id":            sid,
                    "gene_id":              gene_id,
                    "gene_name":            gene_name,
                    "chrom":                chrom,
                    "start":                g_start,
                    "end":                  g_end,
                    "crr":                  crr_val,
                    "segment_cn2_fraction": round(float(elev_frac), 4),
                })

            call_id = _REF_LOOKUP.get((chrom, g_start, g_end))
            if call_id is not None:
                # For reference genes, compute elevated-CN fraction for ALL samples.
                elev_fracs = np.full(n, np.nan)
                for i in range(n):
                    sid = sample_ids[i]
                    elev_fracs[i] = _elevated_cn_fraction(
                        cn_elevated_by_sid.get(sid, []), gene_bin_starts,
                    )

                for i in range(n):
                    crr_val  = float(crr_all[i])   if finite_crr[i]               else None
                    frac_val = float(elev_fracs[i]) if np.isfinite(elev_fracs[i]) else None
                    if not ok[i]:
                        cn = -1
                    elif frac_val is not None and frac_val >= min_gene_coverage_frac:
                        cn = 2
                    else:
                        cn = 1
                    ref_gene_rows.append({
                        "sample_id":            sample_ids[i],
                        "call_id":              call_id,
                        "cn":                   cn,
                        "crr":                  crr_val,
                        "segment_cn2_fraction": round(frac_val, 4) if frac_val is not None else None,
                    })

            n_genes_done += 1
            if n_genes_done % 500 == 0:
                print(f"  {n_genes_done}/{len(genes)} genes | elapsed {time.time() - t0:.0f}s",
                      flush=True)

        del cr_chrom, cn_elevated_by_sid
        gc.collect()

    print(f"Done | {len(genes):,} genes processed | elapsed {time.time() - t0:.0f}s", flush=True)

    del copy_ratios, sanity_ok, segments
    gc.collect()

    # Write genome-wide amplifications (long format, CN>=2 only).
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

    # Write gene_calls.tsv (wide format for 4 reference genes, for evaluation).
    ref_df = pd.DataFrame(ref_gene_rows)
    # Aggregate PM2_PM3 sub-genes: CN=max, CRR=max, fraction=max.
    ref_agg = (
        ref_df.groupby(["sample_id", "call_id"])
              .agg(
                  cn=("cn", "max"),
                  crr=("crr", lambda x: x.max(skipna=True)),
                  segment_cn2_fraction=("segment_cn2_fraction", lambda x: x.max(skipna=True)),
              )
              .reset_index()
    )
    cn_wide   = ref_agg.pivot(index="sample_id", columns="call_id", values="cn")
    crr_wide  = ref_agg.pivot(index="sample_id", columns="call_id", values="crr")
    frac_wide = ref_agg.pivot(index="sample_id", columns="call_id", values="segment_cn2_fraction")
    cn_wide.columns.name  = None
    crr_wide.columns.name = None
    frac_wide.columns.name = None

    ref_ids = list(dict.fromkeys(g["call_id"] for g in REFERENCE_GENES))
    missing = [g for g in ref_ids if g not in cn_wide.columns]
    if missing:
        print(f"WARNING: reference gene(s) {missing} not matched by GFF coordinates — "
              "columns will be NaN in gene_calls.tsv", flush=True)
        for g in missing:
            cn_wide[g]   = pd.NA
            crr_wide[g]  = np.nan
            frac_wide[g] = np.nan

    crr_wide.columns  = [f"crr_{g}"                  for g in crr_wide.columns]
    frac_wide.columns = [f"segment_cn2_fraction_{g}" for g in frac_wide.columns]
    gene_calls_wide = pd.concat([cn_wide, crr_wide, frac_wide], axis=1).reset_index()
    gene_path = os.path.join(out_dir, "gene_calls.tsv")
    gene_calls_wide.to_csv(gene_path, sep="\t", index=False)
    print(f"Saved reference gene calls ({len(gene_calls_wide):,} samples) → {gene_path}", flush=True)
    del ref_gene_rows, ref_df, ref_agg, cn_wide, crr_wide, frac_wide, gene_calls_wide
    gc.collect()

"""
Evaluation — version 05: adds masked CNV diagnostics to v04.

Extends v04 with a MASKED CNV DIAGNOSTICS section that answers:
  "What proportion of ground-truth CNV samples is the model calling as -1
   (no-call / failed), and so excluding from FNR entirely?"

The v04 FNR is computed on the evaluable subset (gt != -1 AND pred != -1).
Samples where the model fails (pred == -1) are silently excluded. If those
failures are non-random — e.g. if weak African amplifications are
disproportionately falling through — then FNR understates total recall loss.

New metrics per gene (overall only; stratified by population/year is left
for future versions if the signal is large enough to warrant it):

  gt_pos_n        — total ground-truth CNV-positive samples
  masked_cnv_n    — gt positives where pred == -1 (excluded from FNR)
  masked_cnv_rate — masked_cnv_n / gt_pos_n
  total_fnr       — (FN_in_eval + masked_cnv_n) / gt_pos_n
                    The inclusive miss rate; if equal to FNR the -1 calls are
                    not biased toward CNV samples. If higher, the -1 calls are
                    silently hiding true misses.

DESIGN NOTE FOR FUTURE EVAL VERSIONS
--------------------------------------
Evaluation modules are the primary tool for detecting blindspots and near-misses.
Each new experiment should consider:
  1. Are there near-miss FNs or FPs that a new diagnostic would expose?
  2. Are existing diagnostics still relevant, or has the calling architecture
     changed enough to make them misleading?
  3. Every metric added must have a clear diagnostic interpretation tied to a
     specific hypothesis about what could go wrong. Scientific defensibility
     is the bar — if you cannot explain what a metric tells you and how you
     would act on it, do not add it.

Reads from cfg:
    pf9_gt_path                  — path to the Pf8/Pf9 curated-coverage TSV (required)
    pf9_meta_path                — path to sample metadata TSV with Population and Year
                                   columns (optional; stratified tables omitted if absent)
    eval_min_group_n             — minimum evaluable samples for a (Year, Population)
                                   group to be printed; default 10
    cnv_min_gene_coverage_fraction — calling threshold; used in near-miss section to
                                   compute recovery counts (default 0.50)
"""

import datetime
import gc
import os

import pandas as pd


GENES    = ["CRT", "GCH1", "MDR1", "PM2_PM3"]
GT_COL   = "{gene}_curated_coverage_only"
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
Q_LABELS  = ["p10", "p25", "p50", "p75", "p90"]


def run_evaluation(out_dir, cfg):
    """Evaluate gene calls against the Pf8 GATK ground truth, plus segment,
    near-miss, and masked-CNV diagnostics.

    Config keys read
    ----------------
    pf9_gt_path                    — path to ground truth TSV (required)
    pf9_meta_path                  — path to sample metadata TSV (optional)
    eval_min_group_n               — min n_eval for (Year, Population) groups (default 10)
    cnv_min_gene_coverage_fraction — calling threshold used for near-miss counts (default 0.50)

    Loads
    -----
    out_dir/gene_calls.tsv    — wide-format gene calls from run_cnv_calls
    out_dir/segments.parquet  — HMM segments for segment diagnostics

    Writes
    ------
    out_dir/evaluation.txt — human-readable metrics report

    Returns
    -------
    dict {"genes": gene_results, "crr": crr_results}
    """
    from sklearn.metrics import matthews_corrcoef  # noqa: PLC0415

    pf9_gt_path    = cfg["pf9_gt_path"]
    pf9_meta_path  = cfg.get("pf9_meta_path")   # optional — None if absent
    min_group_n    = int(cfg.get("eval_min_group_n", 10))
    call_threshold = float(cfg.get("cnv_min_gene_coverage_fraction", 0.50))

    def _cn_to_gt(cn):
        return cn.map(lambda v: -1 if v == -1 else (1 if v > 1 else 0))

    def _metrics(gt, pred_gt):
        """Compute classification metrics for one (gene, group) slice."""
        n         = len(gt)
        eval_mask = (gt != -1) & (pred_gt != -1)
        n_eval    = int(eval_mask.sum())

        # Masked CNV: ground-truth positives that the model calls as -1.
        # These are excluded from FNR, so they are invisible to the standard metric.
        gt_pos_n      = int((gt == 1).sum())
        masked_cnv_n  = int(((gt == 1) & (pred_gt == -1)).sum())

        m = {
            "n":                n,
            "pf9_missing_rate": round((gt == -1).sum() / n, 2) if n else None,
            "new_missing_rate": round((pred_gt == -1).sum() / n, 2) if n else None,
            # delta: model failures on samples PF9 considers callable
            "delta":            round(((gt != -1) & (pred_gt == -1)).sum() / n, 2) if n else None,
            "n_eval":           n_eval,
            "mcc": None, "fnr": None, "ppv": None,
            "gt_pos_n":         gt_pos_n,
            "masked_cnv_n":     masked_cnv_n,
            "masked_cnv_rate":  round(masked_cnv_n / gt_pos_n, 2) if gt_pos_n > 0 else None,
            "total_fnr":        None,  # filled below once fn_count is known
        }
        fn_count = 0
        if n_eval > 0:
            gt_e, pr_e = gt[eval_mask], pred_gt[eval_mask]
            y_true = (gt_e > 0).astype(int)
            y_pred = (pr_e > 0).astype(int)
            if y_true.nunique() == 2:
                m["mcc"] = round(float(matthews_corrcoef(y_true, y_pred)), 2)
            amp_mask = gt_e == 1
            tp = ((gt_e == 1) & (pr_e == 1)).sum()
            fp = ((gt_e == 0) & (pr_e == 1)).sum()
            if amp_mask.sum() > 0:
                fn_count   = int((pr_e[amp_mask] == 0).sum())
                m["fnr"]   = round(fn_count / amp_mask.sum(), 2)
            if (tp + fp) > 0:
                m["ppv"]   = round(float(tp) / float(tp + fp), 2)
        if gt_pos_n > 0:
            m["total_fnr"] = round((fn_count + masked_cnv_n) / gt_pos_n, 2)
        return m

    def _crr_quantiles(series):
        """Return {p10, p25, p50, p75, p90} for a CRR series, or Nones if empty."""
        sub = series.dropna()
        if len(sub) == 0:
            return {"n": 0, **{q: None for q in Q_LABELS}}
        qs = sub.quantile(QUANTILES).round(2).tolist()
        return {"n": len(sub), **dict(zip(Q_LABELS, qs))}

    # gene_calls.tsv is wide: sample_id, CRT, GCH1, MDR1, PM2_PM3, crr_CRT, ...,
    # segment_cn2_fraction_CRT, ...
    _wide_raw = pd.read_csv(os.path.join(out_dir, "gene_calls.tsv"), sep="\t")
    _cast = {g: pd.Int64Dtype() for g in GENES if g in _wide_raw.columns}
    _missing_cols = [g for g in GENES if g not in _wide_raw.columns]
    if _missing_cols:
        print(f"WARNING: gene_calls.tsv missing columns {_missing_cols} — skipping in eval",
              flush=True)
    wide = _wide_raw.astype(_cast).rename(columns={"sample_id": "Sample"})
    del _wide_raw

    gt_cols = ["Sample"] + [GT_COL.format(gene=g) for g in GENES]
    pf9_cnv = pd.read_csv(pf9_gt_path, sep="\t", usecols=gt_cols)
    df      = pf9_cnv.merge(wide, on="Sample")
    del wide, pf9_cnv
    gc.collect()

    has_meta = False
    if pf9_meta_path:
        meta = pd.read_csv(
            pf9_meta_path, sep="\t",
            usecols=["Sample", "Population", "Year"],
        )
        df       = meta.merge(df, on="Sample")
        has_meta = True
        del meta
        gc.collect()

    # ── Compute metrics ──────────────────────────────────────────────────────
    gene_results = {}
    crr_results  = {}

    for gene in GENES:
        if gene not in df.columns:
            print(f"  Skipping {gene}: no calls in gene_calls.tsv", flush=True)
            continue
        gt      = df[GT_COL.format(gene=gene)]
        pred_gt = _cn_to_gt(df[gene])
        crr     = df[f"crr_{gene}"]

        gene_r = {"overall": _metrics(gt, pred_gt)}

        if has_meta:
            gene_r["by_population"] = {
                pop: _metrics(grp[GT_COL.format(gene=gene)], _cn_to_gt(grp[gene]))
                for pop, grp in df.groupby("Population")
            }
            yp_results = {}
            for (year, pop), grp in df.groupby(["Year", "Population"]):
                m = _metrics(grp[GT_COL.format(gene=gene)], _cn_to_gt(grp[gene]))
                if m["n_eval"] >= min_group_n:
                    yp_results[(int(year), pop)] = m
            gene_r["by_year_population"] = yp_results

        gene_results[gene] = gene_r

        # CRR by predicted label
        crr_results[gene] = {
            "by_pred": {
                label: _crr_quantiles(crr[pred_gt == val])
                for val, label in [(-1, "failed"), (0, "pred_normal"), (1, "pred_amp")]
            },
        }
        # CRR by call outcome — evaluable samples only
        eval_mask = (gt != -1) & (pred_gt != -1)
        crr_results[gene]["by_outcome"] = {
            label: _crr_quantiles(crr[eval_mask & (gt == tv) & (pred_gt == pv)])
            for (tv, pv), label in [
                ((0, 0), "TN"), ((0, 1), "FP"), ((1, 0), "FN"), ((1, 1), "TP")
            ]
        }
        # CRR for masked CNV samples (gt==1, pred==-1) — new in v05
        crr_results[gene]["masked_cnv_crr"] = _crr_quantiles(
            crr[(gt == 1) & (pred_gt == -1)]
        )

    # ── Near-miss diagnostics ────────────────────────────────────────────────
    # For each FN sample: how much of the gene body had elevated-CN HMM coverage?
    # Actionable recovery counts: n_fn where fraction >= 0.35 or >= 0.25.
    near_miss_results = {}
    for gene in GENES:
        if gene not in df.columns:
            continue
        frac_col = f"segment_cn2_fraction_{gene}"
        if frac_col not in df.columns:
            near_miss_results[gene] = None
            continue

        gt      = df[GT_COL.format(gene=gene)]
        pred_gt = _cn_to_gt(df[gene])
        eval_mask = (gt != -1) & (pred_gt != -1)
        fn_mask   = eval_mask & (gt == 1) & (pred_gt == 0)
        n_fn = int(fn_mask.sum())

        if n_fn == 0:
            near_miss_results[gene] = {"n_fn": 0}
            continue

        # fillna(0.0): NaN fraction means HMM had no segments for that sample →
        # treat as zero coverage (complete miss)
        fracs = df.loc[fn_mask, frac_col].fillna(0.0)
        q_vals = fracs.quantile(QUANTILES).round(3).tolist()

        near_miss_results[gene] = {
            "n_fn":          n_fn,
            "p10":           q_vals[0],
            "p25":           q_vals[1],
            "p50":           q_vals[2],
            "p75":           q_vals[3],
            "p90":           q_vals[4],
            # recovery counts: how many FNs would be rescued by lowering the threshold
            "n_above_0.35":  int((fracs >= 0.35).sum()),
            "n_above_0.25":  int((fracs >= 0.25).sum()),
            # complete misses: HMM assigned zero elevated-CN coverage to this gene
            "n_zero_frac":   int((fracs == 0.0).sum()),
        }

    del df
    gc.collect()

    # Genes actually evaluated
    eval_genes = [g for g in GENES if g in gene_results]

    # ── Segment diagnostics ──────────────────────────────────────────────────
    seg_diag = _segment_diagnostics(out_dir)

    # ── Format text report ───────────────────────────────────────────────────
    W = 64

    def _fmt(v):
        return f"{v:.2f}" if isinstance(v, float) else ("N/A" if v is None else str(v))

    def _crr_row(label, d):
        qs = " ".join(f"{_fmt(d[q]):>5}" for q in Q_LABELS)
        return f"  {label:<12} {d['n']:>6}  {qs}"

    lines = [
        "=" * W,
        "GUIDANCE",
        "-" * W,
        "FNR (primary metric): fraction of true amplifications called as normal.",
        "  Cross-check with CRR BY OUTCOME below:",
        "    FN p50 >> 1.0 → signal present, HMM discarding it.",
        "      Fix: lower self_transition or adjust state initialisation.",
        "    FN p50 ≈ 1.0  → weak signal, likely upstream coverage issue.",
        "PPV: precision on positive calls. Low PPV is acceptable if FNR is",
        "  the priority. Ground truth is imperfect — some apparent FPs may",
        "  be real signal the GATK pipeline missed.",
        "delta: model-added missingness on PF9-callable samples. High delta",
        "  with CRR[failed] p90 > 1.5 suggests rescuable -1 calls.",
        "MCC: N/A means only one class present (common for CRT).",
        "High FNR may indicate the model learned on too many CNV-positive samples,",
        "  causing the CNV read-count profile to be treated as the default (normal)",
        "  state. In that case the model under-flags amplifications because it has",
        "  been trained to expect elevated counts.",
        "Ground truth source: Pf8 GATK-based pipeline",
        "  (https://github.com/malariagen/malariagen-pf8-cnv-calling/tree/master)",
        "  using DetermineContigPloidy, GermlineCNVCaller, and",
        "  PostprocessGermlineCNVCalls. When results diverge from ground truth,",
        "  consider what this established pipeline does that the current model",
        "  struggles to replicate, and where the model may already be doing better.",
        "",
        "SEGMENT DIAGNOSTICS GUIDANCE",
        "-" * W,
        "callability: fraction of samples where HMM ever left CN1.",
        "  (P. falciparum is haploid — CN1 is the normal baseline, not CN2.)",
        "  Low → model too conservative / VAE signal too weak.",
        "  Compare across experiments to track whether architecture changes",
        "  improve the VAE's ability to flag real variation.",
        "transitions (p10–p90): distribution of within-chrom CN state changes per sample.",
        "  Low p50 → HMM almost never segments; high p90 → HMM too jumpy.",
        "  Do NOT over-optimise these numbers — they are diagnostic, not targets.",
        "",
        "NEAR-MISS DIAGNOSTICS GUIDANCE",
        "-" * W,
        "segment_cn2_fraction: fraction of gene bins in CN>=2 HMM state (confidence-",
        "  filtered at cnv_min_confidence). Includes CN=3, CN=4, etc.",
        "  For FN samples:",
        "    frac ≈ 0           → HMM gave no coverage; parameter tuning alone won't fix.",
        "    0 < frac < 0.25    → partial coverage, well below threshold.",
        "    frac in [0.25,0.50)→ near-miss; lowering cnv_min_gene_coverage_fraction",
        "                         to 0.25 or 0.35 would recover these.",
        "  Recovery counts tell you how many FNs a threshold change would rescue.",
        "  Always ask: if we lower the threshold to recover these FNs, how many new",
        "  FPs would we introduce? Check FP segment_cn2_fraction distribution too.",
        "",
        "MASKED CNV DIAGNOSTICS GUIDANCE",
        "-" * W,
        "masked_cnv_rate: proportion of ground-truth CNV samples the model calls as",
        "  -1 (no-call / quality filter failure). These samples are excluded from FNR",
        "  before the metric is computed, so FNR understates total recall loss when",
        "  masked_cnv_rate is nonzero.",
        "total_fnr: (FN_in_eval + masked_cnv_n) / gt_pos_n — the inclusive miss rate.",
        "  If total_fnr ≈ FNR, -1 calls are unbiased (CNV samples fail at the same",
        "  rate as normal samples). If total_fnr >> FNR, -1 calls are disproportionately",
        "  hitting CNV samples — a systematic blindspot beyond what FNR captures.",
        "  Masked CRR: CRR distribution for masked CNV samples. If masked CRR p50 >> 1.0",
        "  the signal is present but quality filters are blocking the call — rescuable",
        "  with threshold changes. If masked CRR p50 ≈ 1.0, coverage is genuinely low.",
        "",
        "=" * W,
        "Experiment evaluation",
        f"Generated : {datetime.datetime.utcnow().isoformat()}",
        f"Out dir   : {out_dir}",
        "Ground truth: Pf8 GATK pipeline (DetermineContigPloidy / GermlineCNVCaller)",
        "  https://github.com/malariagen/malariagen-pf8-cnv-calling/tree/master",
        "=" * W,
        "",
        "OVERALL",
        "-" * W,
        f"{'Gene':<10} {'MCC':>5} {'FNR':>5} {'PPV':>5} {'call_rate':>10} {'n_eval':>8}",
        "-" * W,
    ]
    for gene in eval_genes:
        m = gene_results[gene]["overall"]
        call_rate = round(1.0 - (m["new_missing_rate"] or 0.0), 2)
        lines.append(
            f"{gene:<10} {_fmt(m['mcc']):>5} {_fmt(m['fnr']):>5} "
            f"{_fmt(m['ppv']):>5} {_fmt(call_rate):>10} {m['n_eval']:>8}"
        )

    lines += [
        "", "MISSINGNESS", "-" * W,
        f"{'Gene':<10} {'pf9_miss':>10} {'model_miss':>12} {'delta':>8}", "-" * W,
    ]
    for gene in eval_genes:
        m = gene_results[gene]["overall"]
        lines.append(
            f"{gene:<10} {_fmt(m['pf9_missing_rate']):>10} "
            f"{_fmt(m['new_missing_rate']):>12} {_fmt(m['delta']):>8}"
        )

    # Masked CNV diagnostics section — new in v05
    lines += [
        "",
        "MASKED CNV  (ground-truth positives called -1; excluded from FNR above)",
        "-" * W,
        f"  {'Gene':<10} {'gt_pos':>8} {'masked_n':>10} {'masked_rate':>13}"
        f" {'total_fnr':>10}",
        f"  {'-' * (W - 2)}",
    ]
    for gene in eval_genes:
        m = gene_results[gene]["overall"]
        lines.append(
            f"  {gene:<10} {m['gt_pos_n']:>8} {m['masked_cnv_n']:>10}"
            f" {_fmt(m['masked_cnv_rate']):>13} {_fmt(m['total_fnr']):>10}"
        )
    lines += [
        "  total_fnr = (FN_in_eval + masked_cnv_n) / gt_pos_n.",
        "  If total_fnr >> FNR: -1 calls are disproportionately hitting CNV samples.",
    ]

    # Masked CNV CRR distribution
    lines += [
        "",
        "MASKED CNV CRR  (CRR for gt-positive samples called -1)",
        "-" * W,
        f"  {'Gene':<12} {'n':>6}  {'p10':>5} {'p25':>5} {'p50':>5} {'p75':>5} {'p90':>5}",
    ]
    for gene in eval_genes:
        d = crr_results[gene]["masked_cnv_crr"]
        lines.append(_crr_row(gene, d))
    lines += [
        "  p50 >> 1.0 → signal present, quality filters blocking; rescuable.",
        "  p50 ≈ 1.0  → genuinely low coverage; calling harder for these samples.",
    ]

    # CRR sections
    q_header = "  " + " ".join(f"{q:>5}" for q in Q_LABELS)
    lines += ["", "CRR BY PREDICTED LABEL  (CRR = gene/flank copy ratio)", "-" * W,
              f"  {'label':<12} {'n':>6}  {q_header.strip()}"]
    for gene in eval_genes:
        lines.append(f"  — {gene}")
        for label, d in crr_results[gene]["by_pred"].items():
            lines.append(_crr_row(label, d))

    lines += [
        "",
        "CRR BY CALL OUTCOME  (evaluable samples only)",
        "  FN p50 >> 1.0 = HMM discarding real signal; FN p50 ≈ 1.0 = weak signal",
        "-" * W,
        f"  {'outcome':<12} {'n':>6}  {q_header.strip()}",
    ]
    for gene in eval_genes:
        lines.append(f"  — {gene}")
        for label, d in crr_results[gene]["by_outcome"].items():
            lines.append(_crr_row(label, d))

    # Near-miss diagnostics section
    has_near_miss = any(
        near_miss_results.get(g) is not None and near_miss_results[g].get("n_fn", 0) >= 0
        for g in eval_genes
    )
    if has_near_miss:
        lines += [
            "",
            "NEAR-MISS DIAGNOSTICS  (FN segment_cn2_fraction at calling confidence)",
            f"  calling threshold = {call_threshold:.2f}",
            "-" * W,
            f"  {'Gene':<10} {'n_FN':>6}  {'p10':>5} {'p25':>5} {'p50':>5} {'p75':>5} {'p90':>5}"
            f"  {'>=0.35':>7} {'>=0.25':>7} {'frac=0':>7}",
            f"  {'-' * (W - 2)}",
        ]
        for gene in eval_genes:
            nm = near_miss_results.get(gene)
            if nm is None:
                lines.append(f"  {gene:<10}  (segment_cn2_fraction not in gene_calls.tsv)")
                continue
            n_fn = nm["n_fn"]
            if n_fn == 0:
                lines.append(f"  {gene:<10} {'0':>6}  (no FNs)")
                continue
            def _f3(v):
                return f"{v:.3f}" if v is not None else "  N/A"
            lines.append(
                f"  {gene:<10} {n_fn:>6}  "
                f"{_f3(nm['p10']):>5} {_f3(nm['p25']):>5} {_f3(nm['p50']):>5} "
                f"{_f3(nm['p75']):>5} {_f3(nm['p90']):>5}  "
                f"{nm['n_above_0.35']:>7} {nm['n_above_0.25']:>7} {nm['n_zero_frac']:>7}"
            )
        lines += [
            f"  >=0.35 / >=0.25: FNs recovered by lowering threshold to 0.35 / 0.25.",
            f"  frac=0: HMM gave no elevated-CN coverage — parameter tuning won't fix.",
        ]

    # BY POPULATION (overall, one table per gene)
    if any("by_population" in gene_results[g] for g in eval_genes):
        for gene in eval_genes:
            if "by_population" not in gene_results[gene]:
                continue
            lines += [
                "", f"BY POPULATION — {gene}", "-" * W,
                f"  {'Population':<14} {'MCC':>5} {'FNR':>5} {'PPV':>5} {'n_eval':>8}",
                f"  {'-' * (W - 2)}",
            ]
            for pop, m in sorted(gene_results[gene]["by_population"].items()):
                lines.append(
                    f"  {pop:<14} {_fmt(m['mcc']):>5} {_fmt(m['fnr']):>5} "
                    f"{_fmt(m['ppv']):>5} {m['n_eval']:>8}"
                )

    # BY (YEAR, POPULATION)
    if any(gene_results[g].get("by_year_population") for g in eval_genes):
        for gene in eval_genes:
            yp = gene_results[gene].get("by_year_population", {})
            if not yp:
                continue
            lines += [
                "",
                f"BY YEAR × POPULATION — {gene}  (n_eval ≥ {min_group_n})",
                "-" * W,
                f"  {'Year':<6} {'Population':<14} {'MCC':>5} {'FNR':>5} {'PPV':>5} {'n_eval':>8}",
                f"  {'-' * (W - 2)}",
            ]
            for (year, pop), m in sorted(yp.items()):
                lines.append(
                    f"  {year:<6} {pop:<14} {_fmt(m['mcc']):>5} {_fmt(m['fnr']):>5} "
                    f"{_fmt(m['ppv']):>5} {m['n_eval']:>8}"
                )

    # Segment diagnostics section
    if seg_diag is not None:
        trans = seg_diag["transitions_percentiles"]
        trans_str = "  ".join(
            f"{lbl}={trans[lbl]:.1f}"
            for lbl in ["p10", "p25", "p50", "p75", "p90"]
        )
        lines += [
            "",
            "SEGMENT DIAGNOSTICS",
            "-" * W,
            f"  callability (samples with any non-CN1 segment): "
            f"{seg_diag['callability']:.3f}",
            f"  within-chrom CN transitions per sample (percentiles):",
            f"    {trans_str}",
            f"  n_samples analysed: {seg_diag['n_samples']}",
        ]
    else:
        lines += ["", "SEGMENT DIAGNOSTICS", "-" * W,
                  "  segments.parquet not found — skipped."]

    out_path = os.path.join(out_dir, "evaluation.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Saved evaluation → {out_path}", flush=True)
    return {"genes": gene_results, "crr": crr_results}


def _segment_diagnostics(out_dir):
    """Return callability and CN transition percentiles from segments.parquet.

    Returns None if the file does not exist (e.g. GCH1-only CRR-fallback runs).
    """
    seg_path = os.path.join(out_dir, "segments.parquet")
    if not os.path.exists(seg_path):
        return None

    segs = pd.read_parquet(seg_path, columns=["sample_id", "chrom", "x0", "cn"])

    # callability: did this sample ever leave CN1?
    # P. falciparum is haploid — CN1 is the normal/baseline state, not CN2.
    sample_has_cnv = segs.groupby("sample_id")["cn"].apply(lambda x: (x != 1).any())
    callability    = float(sample_has_cnv.mean())
    n_samples      = len(sample_has_cnv)

    # within-chromosome CN transitions per sample — reported as percentiles
    segs_s = segs.sort_values(["sample_id", "chrom", "x0"])
    segs_s["prev_cn"] = segs_s.groupby(["sample_id", "chrom"])["cn"].shift(1)
    segs_s["is_trans"] = (
        segs_s["prev_cn"].notna() & (segs_s["cn"] != segs_s["prev_cn"])
    )
    trans_per_sample = segs_s.groupby("sample_id")["is_trans"].sum()
    quantiles        = [0.10, 0.25, 0.50, 0.75, 0.90]
    q_labels         = ["p10", "p25", "p50", "p75", "p90"]
    trans_pcts       = dict(zip(
        q_labels,
        trans_per_sample.quantile(quantiles).round(1).tolist(),
    ))

    del segs, segs_s
    gc.collect()

    return {
        "callability":              round(callability, 3),
        "transitions_percentiles":  trans_pcts,
        "n_samples":                n_samples,
    }

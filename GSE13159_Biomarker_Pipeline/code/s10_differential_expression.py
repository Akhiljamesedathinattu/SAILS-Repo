#!/usr/bin/env python3
"""
PIPELINE STEP 17 : Differential Gene Expression Analysis

One patient group versus all the rest, per gene. Welch's t-test, because group
sizes and within-group variances differ substantially across leukaemia subtypes.
P-values are Benjamini-Hochberg adjusted within each comparison.

SCALE WARNING — read before reporting any effect size from this step.

  The GSE13159 series matrix as distributed by GEO is min-max scaled to [0,1]
  PER SAMPLE. It is not log2 expression, despite what step 01's heuristic
  reports. Verify on your own copy with:

      from common import load_matrix
      v = load_matrix("expr_genes").values
      print(v.min(), v.max())          # 0.0 1.0  -> scaled, not log2

  Consequently the difference of group means computed here is NOT a log2 fold
  change. It is a difference of normalised means, bounded by [-1, 1], and in
  practice tiny: on GSE13159 the median absolute difference is ~0.036 and the
  maximum across all gene x group comparisons is ~0.57. A conventional
  "log2FC >= 1" cutoff can never be met and silently returns nothing.

  The column is named `mean_diff`. There is no log2FC column and there cannot
  be one: recovering a genuine fold change would require reprocessing the raw
  CEL files, since per-sample min-max scaling is not invertible from the
  series matrix alone.

  A `log2FC` alias carrying identical values was written by earlier versions of
  this script, purely so older downstream code kept working. It was a
  misleading name and has been removed. Pass --legacy-log2fc-column to restore
  it temporarily while migrating downstream scripts; the flag warns loudly and
  is intended to be deleted, not kept.

  Because the scale is dataset-specific, the effect-size cutoff defaults to a
  PERCENTILE of the observed distribution rather than a fixed number. Use
  --lfc to override with an absolute value if you prefer.

Groups come from step 09. If those were built from clinical metadata, these
p-values are honest. If they were built from clusters, they are optimistically
biased and the report must say so.

Outputs
  results/10_de_results.csv        every gene x every group
  results/10_de_shortlist.csv      ranked markers per group
  results/10_de_summary.csv
  results/10_heatmap_matrix.csv, 10_heatmap_col_anno.csv, 10_heatmap_row_anno.csv
  results/10_volcano_data.csv
  results/10_effect_size_cutoff.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats

from common import RES_DIR, SEED, load_matrix, bh_fdr, ensure_dirs, log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=None,
                    help="absolute |mean difference| cutoff. Omit to derive it "
                         "from --lfc-percentile instead (recommended: the scale "
                         "is dataset-specific, see the SCALE WARNING above).")
    ap.add_argument("--lfc-percentile", type=float, default=99.0,
                    help="when --lfc is not given, use this percentile of the "
                         "observed |mean difference| distribution (default 99).")
    ap.add_argument("--top-per-group", type=int, default=20)
    ap.add_argument("--heatmap-genes", type=int, default=8)
    ap.add_argument("--heatmap-samples", type=int, default=25)
    ap.add_argument("--legacy-log2fc-column", action="store_true",
                    help="DEPRECATED. Also write a `log2FC` column duplicating "
                         "`mean_diff`, for downstream scripts not yet migrated. "
                         "The values are NOT fold changes. Do not report them.")
    args = ap.parse_args()
    ensure_dirs()

    expr = load_matrix("expr_genes")
    pg = pd.read_csv(RES_DIR / "09_patient_groups.csv").set_index("sample")
    pg = pg.reindex(expr.columns)
    groups = pg["patient_group"].astype(str).to_numpy()
    analysable = pg["analysable"].to_numpy().astype(bool)

    usable = sorted(set(groups[analysable]))
    genes = expr.index.to_numpy()
    X = expr.values.astype(np.float64)
    log(f"testing {len(genes)} genes across {len(usable)} patient groups")
    log(f"group source: {pg['group_source'].iloc[0]}")

    if args.legacy_log2fc_column:
        log("WARNING: --legacy-log2fc-column is set. A `log2FC` column will be "
            "written containing normalised mean differences, NOT fold changes. "
            "Migrate downstream scripts to `mean_diff` and drop this flag.")

    # Report the observed data scale so a wrong cutoff is obvious immediately
    # rather than showing up as an unexplained row of zeros.
    vmin, vmax = float(np.nanmin(X)), float(np.nanmax(X))
    scaled = vmax <= 1.01 and vmin >= -0.01
    log(f"expression range: {vmin:.2f} to {vmax:.2f}"
        + ("  (0-1 scaled — differences are NOT log2 fold changes)"
           if scaled else ""))

    # Per-sample scaling is the specific pathology to check for: it makes each
    # value depend on the extremes of its own array, so a gene is not directly
    # comparable across samples. Detect and record it rather than assume.
    col_min, col_max = X.min(axis=0), X.max(axis=0)
    per_sample_scaled = bool(
        np.allclose(col_min, 0.0, atol=1e-6) and np.allclose(col_max, 1.0, atol=1e-6))
    if per_sample_scaled:
        log("scaling detected: EVERY sample spans exactly [0,1] — the matrix is "
            "min-max scaled per sample, not globally. Report this in Methods.")

    # ---------------------------------------------------------- pass 1: test
    frames = []
    for g in usable:
        m = groups == g
        a, b = X[:, m], X[:, ~m]
        t, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
        p = np.nan_to_num(p, nan=1.0)
        diff = a.mean(axis=1) - b.mean(axis=1)
        cols = {
            "group": g, "gene": genes, "n_in_group": int(m.sum()),
            "mean_in_group": a.mean(axis=1), "mean_rest": b.mean(axis=1),
            "mean_diff": diff,
            "t_stat": np.nan_to_num(t),
            "p_value": p, "fdr": bh_fdr(p),
        }
        if args.legacy_log2fc_column:
            cols["log2FC"] = diff        # deprecated alias, identical values
        frames.append(pd.DataFrame(cols))

    de = pd.concat(frames, ignore_index=True)

    # ------------------------------------------------- effect-size threshold
    if args.lfc is not None:
        cutoff = args.lfc
        rule = "explicit --lfc"
        log(f"effect-size cutoff: |mean_diff| >= {cutoff:g}  (given via --lfc)")
    else:
        cutoff = float(np.percentile(de["mean_diff"].abs(), args.lfc_percentile))
        rule = f"{args.lfc_percentile}th percentile"
        log(f"effect-size cutoff: |mean_diff| >= {cutoff:.4f}  "
            f"({args.lfc_percentile:g}th percentile of observed differences)")

    # Single authoritative record of the effect-size rule, so figures and later
    # steps read one value rather than each recomputing or hardcoding it.
    # Written before the early return below, so the file always exists.
    pd.DataFrame({"cutoff": [cutoff], "fdr": [args.fdr], "rule": [rule],
                  "scale": ["min-max [0,1] per sample" if per_sample_scaled
                            else ("min-max [0,1]" if scaled else "unknown")],
                  "effect_size_is_log2_fold_change": [False]}) \
        .to_csv(RES_DIR / "10_effect_size_cutoff.csv", index=False)

    # ------------------------------------------------- pass 2: apply and log
    summary = []
    for g in usable:
        sub = de[de.group == g]
        sig = sub.fdr < args.fdr
        n_up = int((sig & (sub.mean_diff >= cutoff)).sum())
        n_dn = int((sig & (sub.mean_diff <= -cutoff)).sum())
        summary.append({"group": g, "n_samples": int(sub.n_in_group.iloc[0]),
                        "n_up": n_up, "n_down": n_dn})
        log(f"  {str(g)[:32]:<34} n={int(sub.n_in_group.iloc[0]):>5}  "
            f"up={n_up:>5}  down={n_dn:>5}")

    de.to_csv(RES_DIR / "10_de_results.csv", index=False)
    pd.DataFrame(summary).to_csv(RES_DIR / "10_de_summary.csv", index=False)

    up = de[(de.fdr < args.fdr) & (de.mean_diff >= cutoff)].copy()
    if up.empty:
        q = np.percentile(de["mean_diff"].abs(), [95, 99, 99.5, 99.9])
        log("no significant up-regulated genes at this cutoff.")
        log(f"  |mean_diff| percentiles — 95th {q[0]:.4f}, 99th {q[1]:.4f}, "
            f"99.5th {q[2]:.4f}, 99.9th {q[3]:.4f}")
        log("  try a lower --lfc-percentile, or raise --fdr")
        return

    up["rank_score"] = -np.log10(up.p_value.clip(lower=1e-300)) * up.mean_diff
    short = (up.sort_values(["group", "rank_score"], ascending=[True, False])
               .groupby("group").head(args.top_per_group).reset_index(drop=True))
    short.to_csv(RES_DIR / "10_de_shortlist.csv", index=False)
    log(f"shortlist: {len(short)} genes across {short.group.nunique()} groups")

    missing = sorted(set(usable) - set(short.group.unique()))
    if missing:
        log(f"note: {len(missing)} group(s) contributed no markers: "
            + ", ".join(str(m)[:24] for m in missing))

    # ------------------------------------------------------------- heatmap
    hm = short.groupby("group").head(args.heatmap_genes).drop_duplicates("gene")
    rng = np.random.default_rng(SEED)
    cols = []
    for g in usable:
        idx = np.where(groups == g)[0]
        if idx.size == 0:
            continue
        cols.extend(rng.choice(idx, min(args.heatmap_samples, idx.size),
                               replace=False))
    cols = np.array(sorted(cols))

    sub = expr.loc[hm.gene.to_numpy()].iloc[:, cols]
    sd = sub.std(axis=1).replace(0, 1)
    sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0).round(4) \
       .to_csv(RES_DIR / "10_heatmap_matrix.csv")
    pd.DataFrame({"sample": expr.columns[cols], "group": groups[cols],
                  "leukemia_class": pg["diagnosis"].to_numpy()[cols],
                  "cluster": pg["hclust_cluster"].to_numpy()[cols]}) \
        .to_csv(RES_DIR / "10_heatmap_col_anno.csv", index=False)
    pd.DataFrame({"gene": hm.gene.to_numpy(), "marker_of": hm.group.to_numpy()}) \
        .to_csv(RES_DIR / "10_heatmap_row_anno.csv", index=False)

    # ------------------------------------------------------------- volcano
    volcano_cols = ["group", "gene", "mean_diff", "p_value", "fdr"]
    if args.legacy_log2fc_column:
        volcano_cols.insert(3, "log2FC")
    keep = de.fdr < args.fdr
    pd.concat([de[keep], de[~keep].sample(frac=0.1, random_state=SEED)]) \
        [volcano_cols] \
        .to_csv(RES_DIR / "10_volcano_data.csv", index=False)

    log(f"cutoff used: |mean_diff| >= {cutoff:.4f}, FDR < {args.fdr}")
    log("reminder: mean_diff is a normalised mean difference, not a log2 fold "
        "change — see the SCALE WARNING in this file before writing it up.")


if __name__ == "__main__":
    main()
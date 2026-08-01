#!/usr/bin/env python3
"""
PIPELINE STEPS 6-7 : Variance Filtering -> Z-score Normalization

This order is correct and worth stating in the report: variance filtering must
come BEFORE z-scoring, because z-scoring forces every probe to unit variance and
would make a variance-based ranking meaningless afterwards.

Filtering is two-stage:
  expression  keep probes with signal somewhere in the cohort (95th percentile
              above the matrix median), removing probes that are dark everywhere
  variance    keep the top N by median absolute deviation. MAD rather than
              variance, so a probe that is strongly expressed in one rare
              subtype is not passed over in favour of a merely noisy probe.

Outputs
  work/expr_filtered.{parquet,npz}   log2 values, filtered  (for later steps)
  work/expr_zscore.{parquet,npz}     row z-scored           (for clustering/PCA)
  results/04_filtering_summary.csv
  results/04_probe_variability.csv
"""
import argparse

import numpy as np
import pandas as pd

from common import RES_DIR, load_matrix, save_matrix, ensure_dirs, log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-var", type=int, default=5000)
    ap.add_argument("--drop-outliers", action="store_true",
                    help="exclude samples flagged in results/03_qc_outliers.csv")
    args = ap.parse_args()
    ensure_dirs()

    expr = load_matrix("expr_norm")
    n_start = expr.shape[0]

    if args.drop_outliers:
        f = RES_DIR / "03_qc_outliers.csv"
        if f.exists():
            bad = set(pd.read_csv(f)["sample"].astype(str))
            cols = [c for c in expr.columns if c not in bad]
            log(f"dropping {expr.shape[1] - len(cols)} QC outlier samples")
            expr = expr[cols]

    X = expr.values
    p95 = np.percentile(X, 95, axis=1)
    thr = float(np.median(X))
    expressed = p95 > thr
    expr_e = expr[expressed]
    log(f"expression filter (p95 > {thr:.2f}): {expr_e.shape[0]}/{n_start} kept")

    Xe = expr_e.values
    med = np.median(Xe, axis=1, keepdims=True)
    mad = np.median(np.abs(Xe - med), axis=1)
    pd.DataFrame({"probe_id": expr_e.index, "mad": mad}) \
        .sort_values("mad", ascending=False) \
        .to_csv(RES_DIR / "04_probe_variability.csv", index=False)

    n_top = min(args.top_var, Xe.shape[0])
    top = np.sort(np.argsort(mad)[::-1][:n_top])
    expr_f = expr_e.iloc[top]
    log(f"variance filter: top {n_top} probes by MAD retained")
    save_matrix(expr_f, "expr_filtered")

    # STEP 7: z-score each probe across samples
    sd = expr_f.std(axis=1).replace(0, 1)
    z = expr_f.sub(expr_f.mean(axis=1), axis=0).div(sd, axis=0)
    log(f"z-scored: mean {float(z.values.mean()):.2e}, "
        f"sd {float(z.values.std()):.3f} (expect ~0 and ~1)")
    save_matrix(z.astype(np.float32), "expr_zscore")

    pd.DataFrame({
        "stage": ["loaded", "after_expression_filter", "after_variance_filter"],
        "n_probes": [n_start, expr_e.shape[0], expr_f.shape[0]],
        "n_samples": [expr.shape[1], expr.shape[1], expr.shape[1]],
    }).to_csv(RES_DIR / "04_filtering_summary.csv", index=False)


if __name__ == "__main__":
    main()

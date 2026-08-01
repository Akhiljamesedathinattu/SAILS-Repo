#!/usr/bin/env python3
"""
PIPELINE STEPS 4-5 : Quality Assessment -> Missing Value Analysis

Quality assessment here includes the two transforms the raw data requires before
any distance or correlation is meaningful. They are part of QC, not extra steps:

  log2(x+1)            GSE13159 ships on a linear MAS5 scale. Untransformed,
                       variance scales with intensity and every downstream
                       distance is dominated by the brightest probes.
  quantile normalise   The MILE cohort was assembled across many laboratories.
                       Without between-array normalisation, PC1 encodes array
                       brightness rather than biology.

Missing value analysis follows: Affymetrix arrays rarely carry NAs, so this is a
verification rather than an imputation exercise, but the counts are reported and
any NAs are filled with the probe row mean.

Outputs
  work/expr_norm.{parquet,npz}
  results/03_qc_sample_summary.csv
  results/03_qc_outliers.csv
  results/03_missing_values.csv
"""
import argparse

import numpy as np
import pandas as pd

from common import RES_DIR, load_matrix, save_matrix, ensure_dirs, log


def quantile_normalize(X):
    n_feat, n_samp = X.shape
    ref = np.zeros(n_feat, dtype=np.float64)
    for j in range(n_samp):
        ref += np.sort(X[:, j])
    ref = (ref / n_samp).astype(np.float32)
    rank = np.empty(n_feat, dtype=np.int64)
    for j in range(n_samp):
        order = np.argsort(X[:, j], kind="mergesort")
        rank[order] = np.arange(n_feat)
        X[:, j] = ref[rank]
    return X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-qnorm", action="store_true",
                    help="ablation: skip between-array normalisation")
    ap.add_argument("--outlier-z", type=float, default=5.0)
    args = ap.parse_args()
    ensure_dirs()

    expr = load_matrix("expr_raw")

    # ---- STEP 5: missing value analysis -------------------------------------
    na_per_probe = expr.isna().sum(axis=1)
    na_per_sample = expr.isna().sum(axis=0)
    total_na = int(na_per_probe.sum())
    pd.DataFrame({
        "metric": ["total_missing", "pct_missing", "probes_with_any_missing",
                   "samples_with_any_missing", "worst_probe_missing",
                   "worst_sample_missing"],
        "value": [total_na,
                  round(100 * total_na / expr.size, 6),
                  int((na_per_probe > 0).sum()),
                  int((na_per_sample > 0).sum()),
                  int(na_per_probe.max()), int(na_per_sample.max())],
    }).to_csv(RES_DIR / "03_missing_values.csv", index=False)
    log(f"missing values: {total_na} "
        f"({100 * total_na / expr.size:.4f}% of the matrix)")
    if total_na:
        rowmean = expr.mean(axis=1)
        expr = expr.apply(lambda c: c.fillna(rowmean), axis=0)
        log("missing values imputed with probe row mean")

    # ---- STEP 4: quality assessment -----------------------------------------
    keep = ~expr.index.astype(str).str.upper().str.startswith("AFFX")
    log(f"removing {int((~keep).sum())} AFFX control probes")
    expr = expr[keep]

    X = np.ascontiguousarray(expr.values, dtype=np.float32)
    if np.nanmax(X) > 50:
        log("linear scale detected -> log2(x + 1)")
        np.clip(X, 0, None, out=X)
        X = np.log2(X + 1.0, dtype=np.float32)
    else:
        log("data already log-scaled")

    if args.skip_qnorm:
        log("quantile normalisation SKIPPED (ablation)")
    else:
        log("quantile normalising arrays ...")
        X = quantile_normalize(X)

    expr = pd.DataFrame(X, index=expr.index, columns=expr.columns)

    q = np.percentile(X, [0, 25, 50, 75, 100], axis=0)
    qc = pd.DataFrame({"sample": expr.columns, "min": q[0], "q25": q[1],
                       "median": q[2], "q75": q[3], "max": q[4],
                       "mean": X.mean(axis=0), "sd": X.std(axis=0)})

    rng = np.random.default_rng(42)
    sub = rng.choice(X.shape[0], min(4000, X.shape[0]), replace=False)
    C = np.corrcoef(X[sub].T.astype(np.float64))
    np.fill_diagonal(C, np.nan)
    mc = np.nanmean(C, axis=1)
    qc["mean_correlation"] = mc
    z = (mc - mc.mean()) / (mc.std() or 1.0)
    qc["corr_z"] = z
    qc["outlier"] = z < -args.outlier_z
    qc.to_csv(RES_DIR / "03_qc_sample_summary.csv", index=False)
    qc[qc.outlier].to_csv(RES_DIR / "03_qc_outliers.csv", index=False)

    log(f"mean inter-array correlation {mc.mean():.3f} (min {mc.min():.3f})")
    n_out = int(qc.outlier.sum())
    log(f"{n_out} samples flagged as outliers — FLAGGED, NOT REMOVED. "
        "Inspect results/03_qc_outliers.csv and document your decision.")
    log(f"cached -> {save_matrix(expr, 'expr_norm')}")


if __name__ == "__main__":
    main()

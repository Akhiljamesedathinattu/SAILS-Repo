#!/usr/bin/env python3
"""
PIPELINE STEP 11 : PCA Visualization

PCA on the z-scored filtered matrix, used to visualise the clusters assigned in
step 10. Note the ordering consequence, worth one sentence in the report: because
clinical metadata is not integrated until step 15, this projection can only be
coloured by cluster. Step 15 re-exports the same coordinates with diagnosis
attached, giving the more informative figure.

Outputs
  results/06_pca_scores.csv
  results/06_pca_variance.csv
  results/06_pca_loadings.csv
"""
import argparse

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from common import RES_DIR, WORK_DIR, SEED, load_matrix, ensure_dirs, log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--top-loadings", type=int, default=50)
    args = ap.parse_args()
    ensure_dirs()

    z = load_matrix("expr_zscore")
    M = z.values.T.astype(np.float64)          # samples x probes, already z-scored
    n_pcs = int(min(args.n_pcs, min(M.shape) - 1))

    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=SEED)
    scores = pca.fit_transform(M)
    evr = pca.explained_variance_ratio_
    log(f"PC1-PC3 explain {100 * evr[:3].sum():.1f}% of total variance")

    np.savez_compressed(WORK_DIR / "pca_scores.npz",
                        scores=scores.astype(np.float32),
                        samples=np.asarray(z.columns, dtype=str))

    pd.DataFrame({"pc": np.arange(1, n_pcs + 1), "variance_explained": evr,
                  "cumulative": np.cumsum(evr)}) \
        .to_csv(RES_DIR / "06_pca_variance.csv", index=False)

    out = pd.DataFrame({"sample": z.columns})
    for i in range(min(10, n_pcs)):
        out[f"PC{i+1}"] = scores[:, i]
    ca = RES_DIR / "05_cluster_assignments.csv"
    if ca.exists():
        assign = pd.read_csv(ca)
        out = out.merge(assign, on="sample", how="left")
    out.to_csv(RES_DIR / "06_pca_scores.csv", index=False)

    # which probes drive PC1-PC3 — useful for the report's interpretation section
    load_rows = []
    for i in range(min(3, n_pcs)):
        L = pd.Series(pca.components_[i], index=z.index)
        top = L.abs().sort_values(ascending=False).head(args.top_loadings).index
        load_rows.append(pd.DataFrame({"pc": i + 1, "probe_id": top,
                                       "loading": L.loc[top].to_numpy()}))
    pd.concat(load_rows, ignore_index=True).to_csv(
        RES_DIR / "06_pca_loadings.csv", index=False)
    log(f"top loadings -> {RES_DIR/'06_pca_loadings.csv'}")


if __name__ == "__main__":
    main()

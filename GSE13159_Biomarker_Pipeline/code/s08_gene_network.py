#!/usr/bin/env python3
"""
PIPELINE STEP 14 : Gene Correlation Network

Correlation network over the most variable genes, with modules detected by
hierarchical clustering on 1 - |r|.

READ THIS BEFORE INTERPRETING THE OUTPUT. In a cohort spanning many disease
subtypes, gene-gene correlation is dominated by subtype identity: any two genes
that are both high in the same subtype correlate strongly whether or not they are
functionally related. The result is a small number of enormous modules that
simply re-describe the clusters from step 10. The script therefore computes two
versions:

  global    correlation across all samples. Expect large, subtype-driven modules.
  within    correlation computed within each cluster and then averaged, which
            removes the between-subtype signal and leaves co-regulation that
            holds regardless of subtype. This is the more biologically meaningful
            network and usually far sparser.

Report both and say which you are interpreting. If you have time for only one
thing here, WGCNA is the proper tool for this job; this is a defensible
lightweight substitute, not a replacement.

Outputs
  results/08_network_edges_{global,within}.csv
  results/08_network_nodes.csv
  results/08_network_modules.csv
  results/08_network_summary.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from common import RES_DIR, load_matrix, ensure_dirs, log


def corr_matrix(X):
    C = np.corrcoef(X)
    return np.nan_to_num(C, nan=0.0)


def within_cluster_corr(X, labels):
    """Average of per-cluster correlation matrices, weighted by cluster size."""
    acc = np.zeros((X.shape[0], X.shape[0]))
    wsum = 0.0
    for c in np.unique(labels):
        m = labels == c
        if m.sum() < 10:
            continue
        acc += corr_matrix(X[:, m]) * m.sum()
        wsum += m.sum()
    return acc / wsum if wsum else acc


def edges_from(C, genes, thr):
    iu = np.triu_indices_from(C, k=1)
    r = C[iu]
    keep = np.abs(r) >= thr
    return pd.DataFrame({
        "gene_a": genes[iu[0][keep]], "gene_b": genes[iu[1][keep]],
        "r": np.round(r[keep], 4),
        "sign": np.where(r[keep] > 0, "positive", "negative"),
    })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-genes", type=int, default=300,
                    help="most variable genes to include (keeps the network readable)")
    ap.add_argument("--threshold", type=float, default=0.7,
                    help="minimum |r| for an edge")
    ap.add_argument("--n-modules", type=int, default=8)
    args = ap.parse_args()
    ensure_dirs()

    expr = load_matrix("expr_genes")
    X = expr.values.astype(np.float64)
    med = np.median(X, axis=1, keepdims=True)
    mad = np.median(np.abs(X - med), axis=1)
    top = np.argsort(mad)[::-1][:min(args.n_genes, X.shape[0])]
    sub = expr.iloc[np.sort(top)]
    genes = sub.index.to_numpy()
    Xs = sub.values.astype(np.float64)
    log(f"network over the {len(genes)} most variable genes")

    C_glob = corr_matrix(Xs)
    e_glob = edges_from(C_glob, genes, args.threshold)
    e_glob.to_csv(RES_DIR / "08_network_edges_global.csv", index=False)
    log(f"global network: {len(e_glob)} edges at |r| >= {args.threshold} "
        f"(density {2 * len(e_glob) / (len(genes) * (len(genes) - 1)):.3f})")

    ca = RES_DIR / "05_cluster_assignments.csv"
    e_within = pd.DataFrame()
    if ca.exists():
        assign = pd.read_csv(ca).set_index("sample").reindex(sub.columns)
        labels = assign["hclust_cluster"].to_numpy()
        C_within = within_cluster_corr(Xs, labels)
        e_within = edges_from(C_within, genes, args.threshold)
        e_within.to_csv(RES_DIR / "08_network_edges_within.csv", index=False)
        log(f"within-cluster network: {len(e_within)} edges "
            f"({100 * len(e_within) / max(len(e_glob), 1):.1f}% of the global count)")
        log("the drop is the subtype-driven inflation being removed")
    else:
        C_within = C_glob

    # modules on the global network, since that is what the report will show
    D = 1.0 - np.abs(C_glob)
    np.fill_diagonal(D, 0.0)
    Zl = linkage(squareform(D, checks=False), method="average")
    modules = fcluster(Zl, t=args.n_modules, criterion="maxclust")

    deg_g = pd.concat([e_glob.gene_a, e_glob.gene_b]).value_counts()
    deg_w = (pd.concat([e_within.gene_a, e_within.gene_b]).value_counts()
             if len(e_within) else pd.Series(dtype=int))

    nodes = pd.DataFrame({
        "gene": genes, "module": modules,
        "degree_global": deg_g.reindex(genes).fillna(0).astype(int).to_numpy(),
        "degree_within": deg_w.reindex(genes).fillna(0).astype(int).to_numpy(),
        "mad": mad[np.sort(top)],
        "mean_expression": Xs.mean(axis=1),
    }).sort_values("degree_global", ascending=False)
    nodes.to_csv(RES_DIR / "08_network_nodes.csv", index=False)

    mod = (nodes.groupby("module")
                .agg(n_genes=("gene", "size"),
                     mean_degree=("degree_global", "mean"),
                     hub_gene=("gene", "first"))
                .reset_index())
    mod.to_csv(RES_DIR / "08_network_modules.csv", index=False)
    log("modules (global network):")
    for _, r in mod.iterrows():
        log(f"    module {int(r.module)}: {int(r.n_genes):>4} genes, "
            f"hub={r.hub_gene}")

    pd.DataFrame([{
        "n_genes": len(genes), "threshold": args.threshold,
        "edges_global": len(e_glob), "edges_within": len(e_within),
        "n_modules": int(nodes.module.nunique()),
        "largest_module": int(mod.n_genes.max()),
    }]).to_csv(RES_DIR / "08_network_summary.csv", index=False)

    if len(mod) and mod.n_genes.max() > 0.5 * len(genes):
        log("WARNING: one module holds over half the genes. That is the "
            "subtype-domination effect described in the docstring — interpret "
            "the within-cluster network instead.")


if __name__ == "__main__":
    main()

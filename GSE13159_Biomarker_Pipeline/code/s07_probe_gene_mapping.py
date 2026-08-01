#!/usr/bin/env python3
"""
PIPELINE STEPS 12-13 : Probe-to-Gene Mapping -> Representative Gene Identification

Multiple probes interrogate one gene on HG-U133 Plus 2.0. Left uncollapsed, a
single gene appears many times in every downstream shortlist and inflates the
multiple-testing burden. For each gene symbol the probe with the highest mean
expression is kept ("representative probe"), which is the standard max-mean rule.

This runs on the full normalised matrix from step 4, not the variance-filtered
one, so downstream statistics see all expressed genes rather than only the 5000
most variable.

Requires raw/GPL570_annot.csv. Generate it once with:
    Rscript optional_annotate_probes.R

Outputs
  work/expr_genes.{parquet,npz}
  results/07_probe_gene_mapping.csv
  results/07_representative_genes.csv
"""
import argparse

import numpy as np
import pandas as pd

from common import (RES_DIR, load_matrix, save_matrix, load_probe_annotation,
                    ensure_dirs, log)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-expression-filter", action="store_true", default=True)
    args = ap.parse_args()
    ensure_dirs()

    expr = load_matrix("expr_norm")
    ann = load_probe_annotation()
    if ann is None:
        log("ERROR: raw/GPL570_annot.csv not found.")
        log("Run:  Rscript optional_annotate_probes.R")
        log("Without it, GO/KEGG enrichment (steps 18-19) cannot run at all, "
            "because pathway databases are indexed by gene symbol.")
        raise SystemExit(1)

    sym = ann.reindex(expr.index)
    mapped = sym.notna().to_numpy()
    log(f"{int(mapped.sum())}/{len(expr)} probes map to a gene symbol "
        f"({100 * mapped.mean():.1f}%)")

    sub, sym = expr[mapped], sym[mapped]
    means = sub.mean(axis=1).to_numpy()
    order = np.lexsort((-means, sym.to_numpy()))
    sub_s = sub.iloc[order]
    s = sym.to_numpy()[order]
    first = np.concatenate(([True], s[1:] != s[:-1]))

    genes = sub_s[first]
    rep_probe = genes.index.to_numpy()
    genes.index = pd.Index(s[first], name="gene")
    genes = genes.sort_index()
    log(f"collapsed {int(mapped.sum())} probes -> {genes.shape[0]} unique genes")

    n_probes_per_gene = pd.Series(s).value_counts()
    pd.DataFrame({
        "gene": s[first], "representative_probe": rep_probe,
        "n_probes_for_gene": n_probes_per_gene.reindex(s[first]).to_numpy(),
        "mean_expression": means[order][first],
    }).sort_values("n_probes_for_gene", ascending=False) \
      .to_csv(RES_DIR / "07_representative_genes.csv", index=False)

    pd.DataFrame({"probe_id": expr.index,
                  "gene": sym.reindex(expr.index).to_numpy()}) \
        .to_csv(RES_DIR / "07_probe_gene_mapping.csv", index=False)

    if args.min_expression_filter:
        X = genes.values
        p95 = np.percentile(X, 95, axis=1)
        thr = float(np.median(X))
        keep = p95 > thr
        log(f"expression filter (p95 > {thr:.2f}): {int(keep.sum())} genes retained "
            "— this is the enrichment BACKGROUND universe")
        genes = genes[keep]

    log(f"cached -> {save_matrix(genes, 'expr_genes')}")
    log(f"median probes per collapsed gene: {int(n_probes_per_gene.median())}")


if __name__ == "__main__":
    main()

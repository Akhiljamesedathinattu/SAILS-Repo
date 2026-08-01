#!/usr/bin/env python3
"""
PIPELINE STEPS 15-16 : Clinical Metadata Integration -> Patient Group Creation

THE MOST IMPORTANT STEP IN THE PIPELINE, for a reason that is easy to miss.

Patient groups are defined from the CLINICAL METADATA (the curated diagnosis),
not from the clusters found in step 10. Everything downstream - differential
expression, GO/KEGG, biomarkers, the classifier, the ROC curves - is therefore
anchored to external labels rather than to labels this pipeline invented from the
same expression values. That is what keeps the ROC analysis in step 24 from being
circular and meaningless.

Use --group-source cluster only if you deliberately want the unsupervised
variant, and if you do, say plainly in the report that the resulting p-values and
AUCs are optimistically biased because groups and tests share the same data.

This step is also where clusters finally meet diagnoses: the contingency table
and adjusted Rand index here are the validation of steps 8-10.

Outputs
  results/09_sample_metadata.csv
  results/09_patient_groups.csv
  results/09_cluster_vs_diagnosis.csv
  results/09_cluster_agreement.csv
  results/09_pca_with_diagnosis.csv
  results/09_covariate_association.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from common import RES_DIR, WORK_DIR, load_matrix, ensure_dirs, log


def covariate_association(scores, meta, n_pcs=10):
    """Which metadata fields explain each PC? This is the batch check."""
    rows = []
    for col in meta.columns:
        vals = meta[col].astype(str)
        counts = vals.value_counts()
        usable = counts[counts >= 3].index
        if not (2 <= len(usable) <= 60):
            continue
        mask = vals.isin(usable).to_numpy()
        for pc in range(min(n_pcs, scores.shape[1])):
            x, g = scores[mask, pc], vals[mask].to_numpy()
            try:
                h, p = stats.kruskal(*[x[g == u] for u in usable])
            except ValueError:
                continue
            n = len(x)
            rows.append({"covariate": col, "pc": pc + 1, "n_levels": len(usable),
                         "H": h, "p_value": p,
                         "eta_squared": max(0.0, (h - len(usable) + 1) /
                                            (n - len(usable)))})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label-field", default="leukemia_class",
                    help="metadata column defining patient groups")
    ap.add_argument("--group-source", default="diagnosis",
                    choices=["diagnosis", "cluster"],
                    help="diagnosis = external labels (recommended); "
                         "cluster = unsupervised labels (circular)")
    ap.add_argument("--min-group", type=int, default=20)
    args = ap.parse_args()
    ensure_dirs()

    md_file = WORK_DIR / "sample_metadata_raw.csv"
    if not md_file.exists():
        raise FileNotFoundError("run s02_expression_matrix.py first")
    meta = pd.read_csv(md_file, index_col=0)
    expr = load_matrix("expr_genes")
    meta = meta.reindex(expr.columns)
    log(f"metadata integrated for {len(meta)} samples, "
        f"{meta.shape[1]} fields")

    if args.label_field not in meta.columns:
        raise KeyError(f"'{args.label_field}' not in metadata. Available: "
                       f"{list(meta.columns)}")
    meta.to_csv(RES_DIR / "09_sample_metadata.csv")

    assign = pd.read_csv(RES_DIR / "05_cluster_assignments.csv").set_index("sample")
    assign = assign.reindex(expr.columns)

    diagnosis = meta[args.label_field].astype(str)
    vc = diagnosis.value_counts()
    log(f"{vc.size} distinct values in '{args.label_field}'; "
        f"largest={vc.iloc[0]}, smallest={vc.iloc[-1]}")

    if args.group_source == "diagnosis":
        groups = diagnosis.to_numpy()
        note = "external clinical labels"
    else:
        groups = assign["hclust_cluster"].astype(str).to_numpy()
        note = "unsupervised clusters (CIRCULAR - see docstring)"
    log(f"patient groups from: {note}")

    gc = pd.Series(groups).value_counts()
    usable = gc[gc >= args.min_group].index
    analysable = np.isin(groups, usable)
    log(f"{int(analysable.sum())}/{len(groups)} samples in {len(usable)} groups "
        f"with >= {args.min_group} members")

    pd.DataFrame({
        "sample": expr.columns,
        "patient_group": groups,
        "group_source": args.group_source,
        "diagnosis": diagnosis.to_numpy(),
        "hclust_cluster": assign["hclust_cluster"].to_numpy(),
        "kmeans_cluster": assign["kmeans_cluster"].to_numpy(),
        "analysable": analysable,
    }).to_csv(RES_DIR / "09_patient_groups.csv", index=False)

    # ---- validation of steps 8-10 -------------------------------------------
    ct = pd.crosstab(diagnosis, assign["hclust_cluster"])
    ct.to_csv(RES_DIR / "09_cluster_vs_diagnosis.csv")

    truth = diagnosis.astype("category").cat.codes.to_numpy()
    rows = []
    for name in ("hclust_cluster", "kmeans_cluster"):
        lab = assign[name].to_numpy()
        rows.append({"method": name,
                     "k": int(len(np.unique(lab))),
                     "adjusted_rand_index": adjusted_rand_score(truth, lab),
                     "normalized_mutual_info": normalized_mutual_info_score(truth, lab)})
    met = pd.DataFrame(rows)
    met.to_csv(RES_DIR / "09_cluster_agreement.csv", index=False)
    log("\n" + met.to_string(index=False))

    purity = (ct.max(axis=0) / ct.sum(axis=0)).rename("purity")
    dominant = ct.idxmax(axis=0).rename("dominant_diagnosis")
    pd.concat([ct.sum(axis=0).rename("n"), dominant, purity], axis=1) \
        .to_csv(RES_DIR / "09_cluster_composition.csv")
    log(f"mean cluster purity: {purity.mean():.3f}")

    # ---- re-export PCA with diagnosis attached (the informative figure) -----
    npz = WORK_DIR / "pca_scores.npz"
    if npz.exists():
        z = np.load(npz, allow_pickle=False)
        scores, samples = z["scores"], z["samples"].astype(str)
        out = pd.DataFrame({"sample": samples})
        for i in range(min(10, scores.shape[1])):
            out[f"PC{i+1}"] = scores[:, i]
        out = out.merge(pd.read_csv(RES_DIR / "09_patient_groups.csv"),
                        on="sample", how="left")
        out.to_csv(RES_DIR / "09_pca_with_diagnosis.csv", index=False)

        assoc = covariate_association(scores, meta.loc[samples])
        if not assoc.empty:
            assoc.sort_values(["pc", "eta_squared"], ascending=[True, False]) \
                 .to_csv(RES_DIR / "09_covariate_association.csv", index=False)
            top1 = assoc[assoc.pc == 1].head(3)
            log("strongest associations with PC1:")
            for _, r in top1.iterrows():
                log(f"    {r.covariate:<28} eta2={r.eta_squared:.3f}  "
                    f"p={r.p_value:.2e}")
            log("if a technical field rivals diagnosis here, that is a batch "
                "effect and must be reported")


if __name__ == "__main__":
    main()

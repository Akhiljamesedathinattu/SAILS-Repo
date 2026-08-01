#!/usr/bin/env python3
"""
PIPELINE STEPS 23-25 : Consensus Biomarkers -> ROC Analysis -> Biomarker Validation

STEP 23 combines two independent lines of evidence:
    statistical  candidate score from step 20 (differential expression based)
    predictive   SHAP or impurity importance from step 22
A gene strong on only one is weaker than a gene strong on both. They disagree
more often than people expect, and the overlap is the result worth reporting.
Combined with a harmonic mean, so one-sided support cannot be averaged away.

STEP 24 computes a univariate ROC per consensus gene, on the HELD-OUT TEST SPLIT
only, using the same split as step 22. This matters: a ROC computed on the samples
used to select the gene is inflated and tells you nothing.

STEP 25 validates in three tiers, weakest to strongest:
    tier 1  held-out split of the same cohort   (done here)
    tier 2  a multi-gene panel model trained on train, scored on test (done here)
    tier 3  an independent cohort               (s16_external_validation.py)
Only tier 3 is validation in the strict sense. Say so in the report.

Outputs
  results/consensus_biomarkers.csv
  results/biomarker_roc_curves.csv, biomarker_expression.csv
  results/14_panel_validation.csv, 14_validation_tiers.csv
"""
import argparse

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split

from common import RES_DIR, SEED, load_matrix, ensure_dirs, log


def norm_rank(s):
    r = pd.Series(s).rank(method="average", ascending=True)
    return (r / r.max()).to_numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="random_forest",
                    choices=["random_forest", "logistic_regression"])
    ap.add_argument("--top-importance", type=int, default=500)
    ap.add_argument("--top-out", type=int, default=40)
    ap.add_argument("--roc-genes", type=int, default=12)
    ap.add_argument("--panel-size", type=int, default=25)
    ap.add_argument("--test-size", type=float, default=0.3,
                    help="MUST match s13_machine_learning.py")
    args = ap.parse_args()
    ensure_dirs()

    cand = pd.read_csv(RES_DIR / "12_candidate_biomarkers.csv")
    imp = pd.read_csv(RES_DIR / "ml_feature_importance.csv")
    imp = imp[imp.model == args.model].copy()
    if imp.empty:
        raise ValueError(f"no importance rows for '{args.model}'")
    method = imp.method.iloc[0]
    log(f"predictive evidence: {args.model} ({method})")

    imp = imp.sort_values("importance", ascending=False).drop_duplicates("feature")
    imp["importance_rank"] = np.arange(1, len(imp) + 1)
    pred = imp.head(args.top_importance).set_index("feature")

    stat = cand.drop_duplicates("gene").set_index("gene")
    shared = sorted(set(stat.index) & set(pred.index))
    log(f"{len(stat)} statistical candidates, {len(pred)} predictive; "
        f"{len(shared)} carry BOTH")
    if not shared:
        log("no overlap — raise --top-importance, or check that step 09 and "
            "step 13 used the same grouping")
        return

    tbl = pd.DataFrame({
        "gene": shared,
        "target_group": stat.loc[shared, "group"].to_numpy(),
        "log2FC": stat.loc[shared, "log2FC"].to_numpy(),
        "fdr": stat.loc[shared, "fdr"].to_numpy(),
        "specificity": stat.loc[shared, "specificity"].to_numpy(),
        "candidate_score": stat.loc[shared, "candidate_score"].to_numpy(),
        "importance": pred.loc[shared, "importance"].to_numpy(),
        "importance_rank": pred.loc[shared, "importance_rank"].to_numpy(),
        "importance_method": method,
    })
    tbl["de_percentile"] = norm_rank(tbl.candidate_score)
    tbl["importance_percentile"] = norm_rank(tbl.importance)
    tbl["consensus_score"] = (2 * tbl.de_percentile * tbl.importance_percentile /
                              (tbl.de_percentile + tbl.importance_percentile))
    tbl = tbl.sort_values("consensus_score", ascending=False).reset_index(drop=True)
    tbl.insert(0, "rank", np.arange(1, len(tbl) + 1))

    # ---- STEP 24: univariate ROC on the held-out split ----------------------
    expr = load_matrix("expr_genes")
    pg = pd.read_csv(RES_DIR / "09_patient_groups.csv").set_index("sample")
    pg = pg.reindex(expr.columns)
    y = pg["patient_group"].astype(str).to_numpy()
    keep = pg["analysable"].to_numpy().astype(bool)

    idx = np.arange(expr.shape[1])[keep]
    tr_idx, te_idx = train_test_split(idx, test_size=args.test_size,
                                      stratify=y[keep], random_state=SEED)
    log(f"ROC on {len(te_idx)} held-out samples (same split as step 22)")

    aucs, curves, tidy = [], [], []
    grid = np.linspace(0, 1, 200)
    for _, r in tbl.iterrows():
        if r.gene not in expr.index:
            aucs.append(np.nan)
            continue
        x = expr.loc[r.gene].to_numpy()[te_idx]
        yb = (y[te_idx] == str(r.target_group)).astype(int)
        if yb.sum() < 3 or yb.sum() == len(yb):
            aucs.append(np.nan)
            continue
        a = roc_auc_score(yb, x)
        aucs.append(a)
        if len(curves) < args.roc_genes:
            fpr, tpr, _ = roc_curve(yb, x)
            curves.append(pd.DataFrame({"feature": r.gene,
                                        "group": str(r.target_group), "auc": a,
                                        "fpr": grid,
                                        "tpr": np.interp(grid, fpr, tpr)}))
            tidy.append(pd.DataFrame({"feature": r.gene,
                                      "marker_of": str(r.target_group),
                                      "sample": expr.columns[te_idx],
                                      "group": y[te_idx],
                                      "in_group": np.where(yb == 1,
                                                           str(r.target_group), "rest"),
                                      "expression": x}))
    tbl["univariate_auc_test"] = aucs
    tbl.head(args.top_out).to_csv(RES_DIR / "consensus_biomarkers.csv", index=False)
    if curves:
        pd.concat(curves, ignore_index=True).to_csv(
            RES_DIR / "biomarker_roc_curves.csv", index=False)
        pd.concat(tidy, ignore_index=True).to_csv(
            RES_DIR / "biomarker_expression.csv", index=False)

    log("\n" + tbl.head(10)[["rank", "gene", "target_group", "log2FC",
                             "importance_rank", "univariate_auc_test"]]
        .to_string(index=False))
    log(f"median single-gene AUC: {np.nanmedian(tbl.univariate_auc_test):.3f}")

    # ---- STEP 25 tier 2: does the panel alone work? -------------------------
    panel = tbl.gene.head(args.panel_size).tolist()
    panel = [g for g in panel if g in expr.index]
    Xp = expr.loc[panel].values.T
    clf = RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                 class_weight="balanced_subsample",
                                 n_jobs=-1, random_state=SEED)
    clf.fit(Xp[tr_idx], y[tr_idx])
    pr = clf.predict_proba(Xp[te_idx])
    pp = clf.predict(Xp[te_idx])
    try:
        panel_auc = roc_auc_score(y[te_idx], pr, multi_class="ovr",
                                  average="macro", labels=clf.classes_)
    except ValueError:
        panel_auc = np.nan
    panel_bacc = balanced_accuracy_score(y[te_idx], pp)
    pd.DataFrame([{"panel_size": len(panel), "genes": ";".join(panel),
                   "test_macro_auc": panel_auc,
                   "test_balanced_accuracy": panel_bacc}]) \
        .to_csv(RES_DIR / "14_panel_validation.csv", index=False)
    log(f"{len(panel)}-gene panel: test macro AUC {panel_auc:.3f}, "
        f"balanced accuracy {panel_bacc:.3f}")

    full = pd.read_csv(RES_DIR / "ml_metrics.csv")
    fr = full[full.model == args.model]
    tiers = [
        {"tier": 1, "description": "single genes, held-out split",
         "metric": "median univariate AUC",
         "value": float(np.nanmedian(tbl.univariate_auc_test))},
        {"tier": 2, "description": f"{len(panel)}-gene panel, held-out split",
         "metric": "macro AUC", "value": panel_auc},
        {"tier": 2, "description": "all-gene model, held-out split",
         "metric": "macro AUC",
         "value": float(fr.test_macro_auc_ovr.iloc[0]) if not fr.empty else np.nan},
        {"tier": 3, "description": "independent cohort",
         "metric": "macro AUC", "value": np.nan},
    ]
    pd.DataFrame(tiers).to_csv(RES_DIR / "14_validation_tiers.csv", index=False)
    log("tier 3 (independent cohort) is still empty — run "
        "s16_external_validation.py to fill it. Until then, nothing here is "
        "externally validated.")


if __name__ == "__main__":
    main()

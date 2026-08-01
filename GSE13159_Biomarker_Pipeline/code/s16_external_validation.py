#!/usr/bin/env python3
"""
STEP 9 — External validation on an independent cohort (OPTIONAL)

This is the only step that produces a number you can call validation. Everything
before it is internal: clusters, DE genes and model accuracy all come from one
dataset. Applying the frozen model to a cohort it has never seen is the real test.

Usage
    python3 09_validate_external.py --input raw/GSE13164_series_matrix.txt.gz

Suitable cohorts: the MILE stage-II series, or any GPL570 leukaemia series with
diagnoses in its metadata. The class names will not match GSE13159 exactly, so
--map lets you supply a CSV of external_label,training_label pairs.

Honest caveats this script prints for you:
  - a different cohort means a different batch; a performance drop is expected
    and is itself a finding worth reporting
  - genes absent from the external platform are imputed with the TRAINING mean,
    which biases predictions towards the training majority
  - if the external labels use a different vocabulary, unmapped classes are
    excluded rather than guessed at

Outputs
  results/external_metrics.csv
  results/external_predictions.csv
  results/external_confusion.csv
  results/external_roc_curves.csv
"""
import argparse

import numpy as np
import pandas as pd
from sklearn.metrics import (auc, balanced_accuracy_score, confusion_matrix,
                             roc_auc_score, roc_curve)

from common import (RES_DIR, MODEL_DIR, load_matrix, load_probe_annotation,
                    ensure_dirs, log)
from geo import load_series_matrix


def quantile_normalize_to(X, ref):
    """Map each column of X onto a fixed reference distribution by rank."""
    ref = np.sort(np.asarray(ref, dtype=np.float32))
    n_feat = X.shape[0]
    tgt = np.interp(np.linspace(0, 1, n_feat),
                    np.linspace(0, 1, ref.size), ref).astype(np.float32)
    rank = np.empty(n_feat, dtype=np.int64)
    for j in range(X.shape[1]):
        order = np.argsort(X[:, j], kind="mergesort")
        rank[order] = np.arange(n_feat)
        X[:, j] = tgt[rank]
    return X


def collapse(expr, ann):
    sym = ann.reindex(expr.index)
    mask = sym.notna().to_numpy()
    sub, sym = expr[mask], sym[mask]
    means = sub.mean(axis=1).to_numpy()
    order = np.lexsort((-means, sym.to_numpy()))
    sub = sub.iloc[order]
    s = sym.to_numpy()[order]
    first = np.concatenate(([True], s[1:] != s[:-1]))
    out = sub[first]
    out.index = pd.Index(s[first], name="gene")
    return out.sort_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="external series matrix .gz")
    ap.add_argument("--model", default="random_forest")
    ap.add_argument("--label", default="leukemia_class")
    ap.add_argument("--map", default=None,
                    help="CSV with columns external_label,training_label")
    args = ap.parse_args()
    ensure_dirs()

    try:
        import joblib
    except ImportError:
        raise SystemExit("joblib required: pip install joblib")

    model_path = MODEL_DIR / f"{args.model}.joblib"
    if not model_path.exists():
        raise FileNotFoundError(f"{model_path} missing — run step 07 first")
    pipe = joblib.load(model_path)
    classes = list(pipe.classes_)
    log(f"loaded {args.model} trained on {len(classes)} classes")

    train = load_matrix("expr_genes")
    train_genes = train.index.to_numpy()
    train_mean = train.mean(axis=1)
    train_ref = train.iloc[:, 0].to_numpy()

    log(f"reading external cohort {args.input}")
    expr, meta = load_series_matrix(args.input)
    log(f"external: {expr.shape[0]} probes x {expr.shape[1]} samples")

    expr = expr[~expr.index.astype(str).str.upper().str.startswith("AFFX")]
    X = np.ascontiguousarray(expr.values, dtype=np.float32)
    if np.nanmax(X) > 50:
        np.clip(X, 0, None, out=X)
        X = np.log2(X + 1.0, dtype=np.float32)
        log("external data log2-transformed")
    X = quantile_normalize_to(X, train_ref)
    expr = pd.DataFrame(X, index=expr.index, columns=expr.columns)

    ann = load_probe_annotation()
    if ann is not None:
        expr = collapse(expr, ann)
        log(f"collapsed to {expr.shape[0]} genes")

    present = np.isin(train_genes, expr.index.to_numpy())
    log(f"{int(present.sum())}/{len(train_genes)} training genes found "
        f"({100 * present.mean():.1f}%)")
    if present.mean() < 0.5:
        log("WARNING: under half the training genes are present. Predictions "
            "will be dominated by imputed values — interpret with caution.")

    aligned = expr.reindex(train_genes)
    n_imputed = int(aligned.isna().all(axis=1).sum())
    aligned = aligned.T.fillna(train_mean).T
    if n_imputed:
        log(f"{n_imputed} missing genes imputed with the training mean")

    Xe = aligned.values.T
    pred = pipe.predict(Xe)
    proba = pipe.predict_proba(Xe)

    out = pd.DataFrame({"sample": aligned.columns, "predicted": pred,
                        "confidence": proba.max(axis=1)})
    for i, c in enumerate(classes):
        out[f"p_{c}"] = proba[:, i]

    truth = None
    if args.label in meta.columns:
        truth = meta.reindex(aligned.columns)[args.label].astype(str).to_numpy()
        if args.map:
            m = pd.read_csv(args.map)
            mp = dict(zip(m.external_label.astype(str),
                          m.training_label.astype(str)))
            truth = np.array([mp.get(t, t) for t in truth])
        out["true_label"] = truth
        out["mapped_to_training_class"] = np.isin(truth, classes)
    else:
        log(f"no '{args.label}' field in the external metadata — "
            "predictions saved, but no metrics can be computed")

    out.to_csv(RES_DIR / "external_predictions.csv", index=False)
    log(f"predictions -> {RES_DIR/'external_predictions.csv'}")

    if truth is None:
        return

    keep = np.isin(truth, classes)
    log(f"{int(keep.sum())}/{len(truth)} external samples carry a label that "
        "exists in the training set")
    if keep.sum() < 10:
        log("too few comparable samples for meaningful metrics — "
            "consider a --map file to reconcile the label vocabularies")
        return

    yt, yp, pr = truth[keep], pred[keep], proba[keep]
    try:
        macro = roc_auc_score(yt, pr, multi_class="ovr", average="macro",
                              labels=classes)
    except ValueError:
        macro = np.nan

    pd.DataFrame([{
        "cohort": args.input, "model": args.model,
        "n_samples_scored": int(len(truth)),
        "n_samples_comparable": int(keep.sum()),
        "gene_overlap_fraction": float(present.mean()),
        "external_balanced_accuracy": balanced_accuracy_score(yt, yp),
        "external_macro_auc_ovr": macro,
    }]).to_csv(RES_DIR / "external_metrics.csv", index=False)

    labs = sorted(set(yt) | set(yp))
    pd.DataFrame(confusion_matrix(yt, yp, labels=labs), index=labs, columns=labs) \
        .to_csv(RES_DIR / "external_confusion.csv")

    grid = np.linspace(0, 1, 200)
    rows = []
    for i, c in enumerate(classes):
        yb = (yt == c).astype(int)
        if yb.sum() < 3 or yb.sum() == len(yb):
            continue
        fpr, tpr, _ = roc_curve(yb, pr[:, i])
        rows.append(pd.DataFrame({"cohort": "external", "class": c,
                                  "auc": auc(fpr, tpr), "fpr": grid,
                                  "tpr": np.interp(grid, fpr, tpr)}))
    if rows:
        pd.concat(rows, ignore_index=True).to_csv(
            RES_DIR / "external_roc_curves.csv", index=False)

    internal = RES_DIR / "ml_metrics.csv"
    if internal.exists():
        m = pd.read_csv(internal)
        row = m[m.model == args.model]
        if not row.empty:
            log(f"internal test macro AUC {row.test_macro_auc_ovr.iloc[0]:.3f} "
                f"-> external {macro:.3f}")
            log("report the drop, do not hide it: it is the honest measure of "
                "how far these findings generalise")


if __name__ == "__main__":
    main()

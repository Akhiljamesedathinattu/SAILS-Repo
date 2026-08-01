#!/usr/bin/env python3
"""
PIPELINE STEPS 21-22 : Machine Learning -> SHAP Explainability

Predicts the patient groups defined in step 09. If those came from clinical
metadata, this is a legitimate supervised problem and the AUCs in step 24 mean
something.

LEAKAGE CONTROL, which is what makes or breaks this step: feature selection sits
INSIDE a scikit-learn Pipeline. It is therefore refit on the training portion of
every cross-validation fold and never sees the held-out test split. Selecting
genes on the full dataset and then cross-validating is leakage even when the
labels are legitimate, and it is the single most common error in student
biomarker projects. This applies to the per-class selector below too — it is a
proper transformer, refit per fold, not a precomputed gene list.

SHAP is computed for BOTH models: TreeExplainer for the random forest,
LinearExplainer for logistic regression. Both fall back to the model's native
importance (impurity / mean |coefficient|) if SHAP fails, and the real exception
is logged rather than a generic "install shap" hint.

SHAP OUTPUT LAYOUT — the thing that breaks silently between shap versions.

  shap <= 0.4x returned multiclass values as a LIST of per-class arrays, each
  (n_samples, n_features); np.asarray gave (n_classes, n_samples, n_features).
  shap >= 0.5x returns a single ARRAY (n_samples, n_features, n_classes).

  Code that averages "all axes but the last" is correct for the old layout and
  silently produces a per-CLASS vector on the new one — length n_classes instead
  of n_features, which surfaces as an opaque ValueError. _shap_to_nfc() below
  normalises either layout by matching axis lengths, so this cannot drift again.

FEATURE SELECTION — why --select per_class exists.

  The default SelectKBest(f_classif) scores each gene by how well it separates
  ALL classes simultaneously. On GSE13159 that discarded genes the step-12
  specificity analysis had flagged as its best candidates — BCL11B, GATA3,
  PRKCQ, TOX and CTSW among them, all canonical T-cell genes, i.e. exactly the
  markers wanted for T-ALL. A gene that cleanly separates one subtype of 174
  samples from 1909 others earns a mediocre 17-class F score, because most of
  the between-class variance it is scored against is variance it says nothing
  about.

  --select per_class instead runs a one-vs-rest F test for every class and takes
  the union of each class's top genes. Small subtypes then get representation
  they cannot win in a global ranking. Compare the two: if per_class raises the
  step-12 overlap it supports the selection-bias explanation; if it does not,
  the two methods genuinely disagree and that belongs in the write-up.

  RECORDED HONESTLY: an earlier version of this file asserted that averaging
  SHAP across one-vs-rest classifiers was what suppressed the step-12
  candidates. That was tested per class and it was NOT the cause — on
  GSE13159 only 2 of 24 surviving candidates reached the per-class top 50. The
  overlap is genuinely low. The likely reason is redundancy: with thousands of
  correlated genes a marker can be individually informative yet contribute
  nothing once a stronger correlated neighbour is in the model. Univariate tests
  and multivariate attribution answer different questions, and reporting that
  disagreement is more useful than tuning until it disappears.

Outputs
  results/ml_metrics.csv
  results/ml_roc_curves.csv
  results/ml_confusion_<model>.csv
  results/ml_feature_importance.csv      global ranking, one row per gene
  results/ml_shap_per_class.csv          per gene x class, with rank_in_class
  results/ml_selected_features.csv       which genes survived selection
  models/<model>.joblib
"""
import argparse
import warnings

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.feature_selection._base import SelectorMixin
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (auc, balanced_accuracy_score, confusion_matrix,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from common import RES_DIR, MODEL_DIR, SEED, load_matrix, ensure_dirs, log

warnings.filterwarnings("ignore", category=UserWarning)


class PerClassSelectKBest(SelectorMixin, BaseEstimator):
    """Union of the top-scoring genes from a one-vs-rest F test per class.

    Guarantees every class contributes candidate features, which a single
    global F ranking does not: rare-subtype markers lose to genes that split
    the largest groups. The union size is at most k_per_class * n_classes and
    usually much less, since classes share informative genes.
    """

    def __init__(self, k_per_class=150):
        self.k_per_class = k_per_class

    def fit(self, X, y):
        X = np.asarray(X)
        y = np.asarray(y)
        self.n_features_in_ = X.shape[1]
        classes = np.unique(y)
        mask = np.zeros(X.shape[1], dtype=bool)
        k = int(min(self.k_per_class, X.shape[1]))
        for c in classes:
            yb = (y == c).astype(int)
            with np.errstate(invalid="ignore", divide="ignore"):
                s, _ = f_classif(X, yb)
            s = np.nan_to_num(s, nan=0.0, posinf=0.0, neginf=0.0)
            top = np.argpartition(-s, k - 1)[:k]
            mask[top] = True
        self.support_mask_ = mask
        return self

    def _get_support_mask(self):
        return self.support_mask_


def build_models(k, seed, select="global", k_per_class=150):
    def selector():
        if select == "per_class":
            return PerClassSelectKBest(k_per_class=k_per_class)
        return SelectKBest(f_classif, k=k)

    return {
        "logistic_regression": Pipeline([
            ("sel", selector()),
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced",
                                       random_state=seed))]),
        "random_forest": Pipeline([
            ("sel", selector()),
            ("clf", RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                           class_weight="balanced_subsample",
                                           n_jobs=-1, random_state=seed))]),
    }


def roc_table(y_true, proba, classes, model, n=200):
    grid = np.linspace(0, 1, n)
    rows = []
    for i, c in enumerate(classes):
        yb = (y_true == c).astype(int)
        if yb.sum() == 0 or yb.sum() == len(yb):
            continue
        fpr, tpr, _ = roc_curve(yb, proba[:, i])
        rows.append(pd.DataFrame({"model": model, "class": str(c),
                                  "auc": auc(fpr, tpr), "fpr": grid,
                                  "tpr": np.interp(grid, fpr, tpr)}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _shap_to_nfc(sv, n_samples, n_features, n_classes):
    """Normalise any SHAP return layout to (n_samples, n_features, n_classes).

    Raises with the observed shape rather than letting a mis-oriented array
    through, because a wrong orientation produces plausible-looking numbers
    attached to the wrong genes.
    """
    if isinstance(sv, list):
        sv = np.stack([np.asarray(s) for s in sv], axis=-1)
    sv = np.asarray(sv)

    if sv.ndim == 2:
        return sv[:, :, None]

    if sv.ndim == 3:
        if sv.shape == (n_samples, n_features, n_classes):
            return sv
        if sv.shape == (n_classes, n_samples, n_features):
            return np.transpose(sv, (1, 2, 0))
        if sv.shape == (n_samples, n_classes, n_features):
            return np.transpose(sv, (0, 2, 1))

    raise ValueError(
        f"unrecognised SHAP shape {sv.shape}; expected some permutation of "
        f"(n_samples={n_samples}, n_features={n_features}, n_classes={n_classes})")


def shap_importance(name, clf, Xs, kept, classes, background=None):
    """Return (global_df, per_class_df) from SHAP values, or raise."""
    import shap

    if hasattr(clf, "estimators_"):
        sv = shap.TreeExplainer(clf).shap_values(Xs)
    else:
        bg = background if background is not None else Xs
        sv = shap.LinearExplainer(clf, bg).shap_values(Xs)

    sv = _shap_to_nfc(sv, Xs.shape[0], len(kept), len(classes))
    a = np.abs(sv)

    gdf = pd.DataFrame({"model": name, "feature": kept,
                        "importance": a.mean(axis=(0, 2)),
                        "method": "shap_mean_abs"})

    per_class = a.mean(axis=0)                      # (features, classes)
    cls = [str(c) for c in classes][:per_class.shape[1]]
    pdf = (pd.DataFrame(per_class, index=kept, columns=cls)
             .rename_axis("feature").reset_index()
             .melt(id_vars="feature", var_name="class", value_name="mean_abs_shap"))
    pdf.insert(0, "model", name)
    # rank within class so downstream comparisons need not recompute it
    pdf["rank_in_class"] = (pdf.groupby("class")["mean_abs_shap"]
                               .rank(ascending=False, method="min").astype(int))
    return gdf, pdf


def importance_table(name, pipe, feature_names, X_test, use_shap, n_explain=300):
    kept = np.asarray(feature_names)[pipe.named_steps["sel"].get_support()]
    clf = pipe.named_steps["clf"]
    classes = getattr(clf, "classes_", np.array(["0"]))
    per_class = None

    if use_shap:
        try:
            n = min(n_explain, len(X_test))
            Xs = pipe[:-1].transform(X_test[:n])
            bg = Xs[:min(100, len(Xs))]
            gdf, per_class = shap_importance(name, clf, Xs, kept, classes, bg)
            log(f"  SHAP computed for {name} on {n} test samples "
                f"({len(kept)} features x {len(classes)} classes)")
            return gdf, per_class
        except Exception as e:
            log(f"  SHAP failed for {name}: {type(e).__name__}: {e}")
            log("  falling back to the model's native importance")

    if hasattr(clf, "feature_importances_"):
        gdf = pd.DataFrame({"model": name, "feature": kept,
                            "importance": clf.feature_importances_,
                            "method": "impurity"})
    else:
        gdf = pd.DataFrame({"model": name, "feature": kept,
                            "importance": np.abs(clf.coef_).mean(axis=0),
                            "method": "mean_abs_coefficient"})
    return gdf, per_class


def report_candidate_survival(kept):
    """Say plainly which step-12 candidates never reached the model.

    Without this the overlap statistics downstream are unfair: a gene dropped by
    feature selection cannot appear in any importance ranking, so counting it as
    "not confirmed" conflates two different things.
    """
    f = RES_DIR / "12_candidate_biomarkers.csv"
    if not f.exists():
        return
    cand = pd.read_csv(f)
    if "gene" not in cand.columns:
        return
    genes = set(cand["gene"].astype(str))
    kept_set = set(map(str, kept))
    dropped = sorted(genes - kept_set)
    log(f"step-12 candidates: {len(genes)}, surviving selection "
        f"{len(genes & kept_set)}, dropped {len(dropped)}")
    if dropped:
        log("  dropped: " + ", ".join(dropped[:15])
            + (" ..." if len(dropped) > 15 else ""))
        log("  these cannot appear in any importance ranking — exclude them from "
            "overlap denominators or state the denominator explicitly")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-size", type=float, default=0.3)
    ap.add_argument("--k-features", type=int, default=2000)
    ap.add_argument("--select", choices=["global", "per_class"], default="global",
                    help="global = one 17-class F ranking (original behaviour); "
                         "per_class = union of per-class one-vs-rest top genes")
    ap.add_argument("--k-per-class", type=int, default=150,
                    help="genes kept per class when --select per_class")
    ap.add_argument("--cv-folds", type=int, default=5)
    ap.add_argument("--shap", action="store_true")
    ap.add_argument("--shap-samples", type=int, default=300,
                    help="test samples explained. TreeExplainer cost scales "
                         "with samples x features x classes.")
    args = ap.parse_args()
    ensure_dirs()

    expr = load_matrix("expr_genes")
    pg = pd.read_csv(RES_DIR / "09_patient_groups.csv").set_index("sample")
    pg = pg.reindex(expr.columns)
    y_all = pg["patient_group"].astype(str).to_numpy()
    keep = pg["analysable"].to_numpy().astype(bool)

    src = pg["group_source"].iloc[0]
    log(f"labels = patient groups from {src}")
    if src != "diagnosis":
        log("WARNING: labels came from clusters derived from these same features. "
            "The AUCs below are optimistically biased and must be reported as such.")

    X, y = expr.values.T[keep], y_all[keep]
    genes = expr.index.to_numpy()
    log(f"{X.shape[0]} samples, {len(set(y))} groups, {X.shape[1]} genes")

    X_tr, X_te, y_tr, y_te = train_test_split(
        X, y, test_size=args.test_size, stratify=y, random_state=SEED)
    log(f"train {len(y_tr)} / test {len(y_te)} (stratified, seed {SEED})")

    k = min(args.k_features, X.shape[1])
    if args.select == "per_class":
        log(f"feature selection: per-class union, {args.k_per_class} genes x "
            f"{len(set(y))} classes")
    else:
        log(f"feature selection: global {len(set(y))}-class F test, top {k}")

    cv = StratifiedKFold(n_splits=args.cv_folds, shuffle=True, random_state=SEED)
    metrics, rocs, imps, per_cls, sel_rows = [], [], [], [], []

    models = build_models(k, SEED, args.select, args.k_per_class)
    for name, pipe in models.items():
        log(f"fitting {name} ...")
        cvs = cross_val_score(pipe, X_tr, y_tr, cv=cv,
                              scoring="balanced_accuracy", n_jobs=1)
        pipe.fit(X_tr, y_tr)
        pred, proba = pipe.predict(X_te), pipe.predict_proba(X_te)
        classes = pipe.classes_
        try:
            macro = roc_auc_score(y_te, proba, multi_class="ovr",
                                  average="macro", labels=classes)
        except ValueError:
            macro = np.nan

        kept = genes[pipe.named_steps["sel"].get_support()]
        metrics.append({"model": name,
                        "cv_balanced_accuracy_mean": cvs.mean(),
                        "cv_balanced_accuracy_sd": cvs.std(),
                        "test_balanced_accuracy": balanced_accuracy_score(y_te, pred),
                        "test_macro_auc_ovr": macro,
                        "n_train": len(y_tr), "n_test": len(y_te),
                        "n_classes": len(classes),
                        "k_features": len(kept),
                        "selection": args.select,
                        "label_source": src})
        log(f"  CV {cvs.mean():.3f}+/-{cvs.std():.3f} | "
            f"test {balanced_accuracy_score(y_te, pred):.3f} | "
            f"macro AUC {macro:.3f} | {len(kept)} features used")

        sel_rows.append(pd.DataFrame({"model": name, "feature": kept,
                                      "selection": args.select}))

        rt = roc_table(y_te, proba, classes, name)
        if not rt.empty:
            rocs.append(rt)

        gdf, pdf = importance_table(name, pipe, genes, X_te, args.shap,
                                    args.shap_samples)
        imps.append(gdf)
        if pdf is not None:
            per_cls.append(pdf)

        pd.DataFrame(confusion_matrix(y_te, pred, labels=classes),
                     index=classes, columns=classes) \
            .to_csv(RES_DIR / f"ml_confusion_{name}.csv")
        try:
            import joblib
            joblib.dump(pipe, MODEL_DIR / f"{name}.joblib")
        except Exception as e:
            log(f"  model not pickled: {e}")

    report_candidate_survival(
        genes[models["logistic_regression"].named_steps["sel"].get_support()])

    pd.DataFrame(metrics).to_csv(RES_DIR / "ml_metrics.csv", index=False)
    pd.concat(sel_rows, ignore_index=True) \
      .to_csv(RES_DIR / "ml_selected_features.csv", index=False)
    if rocs:
        pd.concat(rocs, ignore_index=True).to_csv(
            RES_DIR / "ml_roc_curves.csv", index=False)
    pd.concat(imps, ignore_index=True) \
      .sort_values(["model", "importance"], ascending=[True, False]) \
      .to_csv(RES_DIR / "ml_feature_importance.csv", index=False)

    if per_cls:
        allc = pd.concat(per_cls, ignore_index=True)
        allc.sort_values(["model", "class", "mean_abs_shap"],
                         ascending=[True, True, False]) \
            .to_csv(RES_DIR / "ml_shap_per_class.csv", index=False)
        log(f"per-class SHAP -> results/ml_shap_per_class.csv "
            f"({len(allc)} rows, includes rank_in_class)")
    else:
        log("no per-class SHAP written (SHAP disabled or failed for both models)")
    log("done")


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
PIPELINE STEP 25 : Final Biomarker Report

Assembles a markdown report from every results CSV, with your real numbers
already filled in. Sections needing human judgement are marked TODO — the
biology, the interpretation and the limitations are yours to write. Everything
countable is filled automatically so you are not transcribing numbers by hand.

Outputs
  results/FINAL_BIOMARKER_REPORT.md
"""
import argparse
from datetime import date

import numpy as np
import pandas as pd

from common import RES_DIR, FIG_DIR, ensure_dirs, log


def rd(name):
    p = RES_DIR / name
    if not p.exists():
        return None
    try:
        df = pd.read_csv(p)
        return None if df.empty else df
    except Exception:
        return None


def md_table(df, cols=None, n=15, floatfmt=3):
    if df is None or df.empty:
        return "_not available — the relevant step has not been run_\n"
    d = df[cols] if cols else df
    d = d.head(n).copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else
                            (f"{v:.2e}" if abs(v) < 1e-3 and v != 0
                             else f"{v:.{floatfmt}f}"))
    head = "| " + " | ".join(str(c) for c in d.columns) + " |"
    sep = "|" + "|".join(["---"] * len(d.columns)) + "|"
    rows = ["| " + " | ".join(str(v) for v in r) + " |"
            for r in d.itertuples(index=False)]
    return "\n".join([head, sep] + rows) + "\n"


def kv(df, key_col, val_col):
    return {} if df is None else dict(zip(df[key_col].astype(str), df[val_col]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--title", default="Candidate Biomarkers in Acute Leukaemia")
    ap.add_argument("--author", default="[your name]")
    args = ap.parse_args()
    ensure_dirs()

    ov = kv(rd("01_dataset_overview.csv"), "property", "value")
    filt = rd("04_filtering_summary.csv")
    csum = rd("05_clustering_summary.csv")
    agree = rd("09_cluster_agreement.csv")
    comp = rd("09_cluster_composition.csv")
    assoc = rd("09_covariate_association.csv")
    desum = rd("10_de_summary.csv")
    net = rd("08_network_summary.csv")
    ensum = rd("enrichment_summary.csv")
    entop = rd("enrichment_top.csv")
    cand = rd("12_candidate_biomarkers.csv")
    mlm = rd("ml_metrics.csv")
    cons = rd("consensus_biomarkers.csv")
    panel = rd("14_panel_validation.csv")
    tiers = rd("14_validation_tiers.csv")
    ext = rd("external_metrics.csv")

    L = []
    A = L.append
    A(f"# {args.title}\n")
    A(f"**Dataset** GSE13159 (MILE study) · **Author** {args.author} · "
      f"**Generated** {date.today().isoformat()}\n")
    A("> Auto-generated from pipeline outputs. Numbers are final; every **TODO** "
      "marks prose you must write yourself.\n")

    A("\n## 1. Dataset\n")
    A(f"- Platform: {ov.get('platform', 'n/a')}")
    A(f"- Samples: {ov.get('n_samples', 'n/a')} · Probes: {ov.get('n_probes', 'n/a')}")
    A(f"- Value scale as downloaded: {ov.get('value_scale', 'n/a')}")
    if filt is not None:
        A(f"- After filtering: {int(filt.n_probes.iloc[-1])} probes retained "
          f"from {int(filt.n_probes.iloc[0])}\n")
    A("\n**TODO** one paragraph on why this cohort suits the question.\n")

    A("\n## 2. Preprocessing and quality control\n")
    mv = kv(rd("03_missing_values.csv"), "metric", "value")
    A(f"- Missing values: {mv.get('total_missing', 'n/a')} "
      f"({mv.get('pct_missing', 'n/a')}% of the matrix)")
    qc = rd("03_qc_sample_summary.csv")
    if qc is not None:
        A(f"- Mean inter-array correlation after normalisation: "
          f"{qc.mean_correlation.mean():.3f}")
        A(f"- Samples flagged as QC outliers: {int(qc.outlier.sum())} "
          "(flagged, not removed)")
    A("- log2 transform and quantile normalisation applied within quality "
      "assessment, because the data ship on a linear scale and the cohort spans "
      "multiple laboratories.\n")

    A("\n## 3. Clustering\n")
    if csum is not None:
        r = csum.iloc[0]
        A(f"- k = {int(r.k_final)} chosen by {r.k_rule}; lowest PAC would "
          f"suggest k = {int(r.k_min_pac)}")
        A(f"- Silhouette {r.silhouette:.3f} · mean consensus "
          f"{r.mean_consensus:.3f}\n")
    A("\nAgreement with curated diagnoses (clusters never saw these labels):\n")
    A(md_table(agree))
    if comp is not None:
        A("\nCluster composition:\n")
        A(md_table(comp, n=25))
    A("\n**TODO** which clusters are pure, which fuse two diagnoses, which "
      "diagnosis splits across clusters, and whether those merges are "
      "biologically sensible.\n")

    if assoc is not None:
        A("\n### Batch check\n")
        A(md_table(assoc[assoc.pc == 1],
                   ["covariate", "pc", "n_levels", "eta_squared", "p_value"], n=6))
        A("\n**TODO** if a technical field rivals diagnosis on PC1, say so "
          "explicitly and discuss the confound.\n")

    if net is not None:
        A("\n## 4. Gene correlation network\n")
        r = net.iloc[0]
        A(f"- {int(r.n_genes)} genes, |r| >= {r.threshold}")
        A(f"- Global edges: {int(r.edges_global)} · within-cluster edges: "
          f"{int(r.edges_within)}")
        A(f"- {int(r.n_modules)} modules, largest holding "
          f"{int(r.largest_module)} genes\n")
        A("\n**TODO** state which network you interpret and why. In a "
          "multi-subtype cohort the global network is dominated by subtype "
          "identity; the within-cluster network is the more meaningful one.\n")

    A("\n## 5. Differential expression\n")
    A(md_table(desum, n=25))
    A("\n**TODO** comment on the relation between group size and gene counts; "
      "small groups have less power but a homogeneous phenotype can still give "
      "very small p-values.\n")

    if ensum is not None:
        A("\n## 6. GO and KEGG enrichment\n")
        A(md_table(ensum, n=20))
        if entop is not None:
            A("\nTop terms:\n")
            A(md_table(entop, ["library", "group", "label", "test", "effect", "fdr"],
                       n=20))
        A("\n**TODO** lead with terms significant under BOTH tests. State that "
          "the background was the filtered gene universe, not the genome.\n")

    A("\n## 7. Candidate biomarkers\n")
    A(md_table(cand, ["rank", "gene", "group", "log2FC", "fdr", "specificity",
                      "candidate_score"], n=20))

    A("\n## 8. Classification\n")
    if mlm is not None:
        A(md_table(mlm, ["model", "cv_balanced_accuracy_mean",
                         "test_balanced_accuracy", "test_macro_auc_ovr",
                         "label_source"]))
        src = str(mlm.label_source.iloc[0])
        if src == "diagnosis":
            A("\nLabels were the curated clinical diagnoses, not the clusters, so "
              "these figures are not circular. Feature selection sat inside the "
              "cross-validation pipeline and never saw the test split.\n")
        else:
            A("\n**WARNING** labels came from clusters derived from these same "
              "features. These figures are optimistically biased and must be "
              "reported as such.\n")

    A("\n## 9. Consensus biomarkers\n")
    if cons is not None:
        A(f"{len(cons)} genes carry both statistical and predictive support.\n")
        A(md_table(cons, ["rank", "gene", "target_group", "log2FC", "fdr",
                          "importance_rank", "univariate_auc_test"], n=20))
        A(f"\nMedian single-gene AUC: "
          f"{np.nanmedian(cons.univariate_auc_test):.3f}\n")
    A("\n**TODO** for each of the top genes, one or two sentences of literature "
      "context. Genes already known to mark that subtype are your positive "
      "control that the pipeline works; unexpected genes are hypotheses, not "
      "discoveries.\n")

    A("\n## 10. Validation\n")
    if panel is not None:
        r = panel.iloc[0]
        A(f"- {int(r.panel_size)}-gene panel on held-out samples: macro AUC "
          f"{r.test_macro_auc:.3f}, balanced accuracy "
          f"{r.test_balanced_accuracy:.3f}\n")
    if tiers is not None:
        A(md_table(tiers))
    if ext is not None:
        A("\n### Independent cohort\n")
        A(md_table(ext))
    else:
        A("\n**No independent cohort was tested.** All figures above are "
          "internal to GSE13159. This is a limitation, not a failure, but it "
          "must be stated plainly: nothing here is externally validated.\n")

    A("\n## 11. Limitations\n")
    A("- One-vs-rest testing makes the reference group heterogeneous, inflating "
      "significance for genes absent from one large competing subtype.")
    A("- Multiple probes per gene were collapsed by max-mean; a different rule "
      "would give a slightly different gene set.")
    A("- Bulk tissue: differences may reflect changing cell-type proportions "
      "rather than regulation within a cell type.")
    A("- No batch covariate was modelled; centre association is measured but "
      "not corrected.")
    if mlm is not None and str(mlm.label_source.iloc[0]) != "diagnosis":
        A("- Groups were defined from the same expression data used to test "
          "them, so all p-values and AUCs are optimistically biased.")
    A("- **TODO** add any dataset-specific limitations you encountered.\n")

    A("\n## Figures\n")
    figs = sorted(p.name for p in FIG_DIR.glob("*.png")) if FIG_DIR.exists() else []
    if figs:
        for f in figs:
            A(f"- `{f}`")
    else:
        A("_No figures found. Run the R scripts._")
    A("\n## Reproducibility\n")
    A("```bash\nexport SAILS_BASE=" + str(RES_DIR.parent) + "\nbash run_all.sh\n```")
    A("\nRandom seed 42 throughout.\n")

    out = RES_DIR / "FINAL_BIOMARKER_REPORT.md"
    out.write_text("\n".join(L), encoding="utf-8")
    log(f"report -> {out}")
    n_todo = "\n".join(L).count("TODO")
    log(f"{n_todo} TODO markers left for you to write")


if __name__ == "__main__":
    main()

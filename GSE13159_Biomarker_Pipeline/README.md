# GSE13159 Biomarker Discovery Pipeline

Implements the 25-step pipeline exactly as specified, in the specified order.
Python does all computation; R draws all figures. They communicate through CSVs
in `results/`.

## Step → script map

| # | Pipeline step | Script |
|---|---|---|
| 1 | Raw GSE13159 Dataset | `s01_data_understanding.py` |
| 2 | Data Understanding | ″ |
| 3 | Expression Matrix Loading | `s02_expression_matrix.py` |
| 4 | Quality Assessment | `s03_quality_assessment.py` |
| 5 | Missing Value Analysis | ″ |
| 6 | Variance Filtering | `s04_filter_normalize.py` |
| 7 | Z-score Normalization | ″ |
| 8 | Hierarchical Clustering | `s05_hierarchical_clustering.py` |
| 9 | Dendrogram Generation | ″ |
| 10 | Cluster Assignment | ″ |
| 11 | PCA Visualization | `s06_pca.py` |
| 12 | Probe-to-Gene Mapping | `s07_probe_gene_mapping.py` |
| 13 | Representative Gene Identification | ″ |
| 14 | Gene Correlation Network | `s08_gene_network.py` |
| 15 | Clinical Metadata Integration | `s09_metadata_groups.py` |
| 16 | Patient Group Creation | ″ |
| 17 | Differential Gene Expression Analysis | `s10_differential_expression.py` |
| 18 | GO Enrichment Analysis | `s11_go_kegg.py` |
| 19 | KEGG Pathway Analysis | ″ |
| 20 | Candidate Biomarker Identification | `s12_candidate_biomarkers.py` |
| 21 | Machine Learning | `s13_machine_learning.py` |
| 22 | SHAP Explainability | ″ |
| 23 | Consensus Biomarkers | `s14_consensus_roc.py` |
| 24 | ROC Analysis | ″ |
| 25 | Biomarker Validation / Final Report | ″, `s16_external_validation.py`, `s15_final_report.py` |

Nothing was reordered. Two things were folded *inside* existing steps:
log2 + quantile normalisation live in step 4 (Quality Assessment), and the
enrichment background is set to the filtered gene universe from step 12.

---

## The one thing that makes this pipeline sound

Step 15 → 16 is *Clinical Metadata Integration → Patient Group Creation*. Patient
groups are therefore built from the **curated diagnosis**, and everything
downstream — differential expression, GO/KEGG, candidate biomarkers, the
classifier, the ROC curves — is anchored to external labels.

That is what keeps step 24 meaningful. Had patient groups been carved out of the
step-10 clusters instead, the classifier would be predicting labels that this
pipeline invented from the same expression values, and the resulting AUC would be
near-perfect and worthless.

The clusters from steps 8–10 remain an independent result. Step 15 compares them
to diagnosis via adjusted Rand index and a contingency table, which validates the
unsupervised half of the pipeline.

`s09_metadata_groups.py --group-source cluster` switches to the circular variant
if you want it deliberately. Every downstream script then prints a warning, and
the generated report inserts a bias caveat automatically.

---

## Quick start

```bash
export SAILS_BASE=/sails/SAILS-Repo/Gene_Expression_Clustering
mkdir -p $SAILS_BASE/raw
# put GSE13159_series_matrix.txt.gz into $SAILS_BASE/raw/

pip install numpy pandas scipy scikit-learn pyarrow joblib
pip install shap                          # optional; step 22 falls back without it

Rscript code/optional_annotate_probes.R   # REQUIRED before step 12
python3 code/fetch_genesets.py            # GO + KEGG libraries for steps 18-19
bash code/run_all.sh
```

### Test it first — one minute, known answers

```bash
export SAILS_BASE=/tmp/bm_test
python3 code/make_synthetic_data.py
bash code/run_all.sh
```

Synthetic data with five classes and planted marker blocks. Verified results on
this test set:

| Check | Expected | Observed |
|---|---|---|
| Cluster vs diagnosis ARI | ≈ 1.0 | 0.986 (Ward), 1.000 (k-means) |
| Mean cluster purity | ≈ 1.0 | 0.995 |
| k chosen by silhouette / PAC | 5 / 5 | 5 / 5 |
| Consensus biomarkers | planted markers | all top 10 correct |
| Median single-gene AUC | high | 0.933 |
| Decoy pathways enriched | 0 | 0 of 100 rank tests |
| Network self-diagnosis | warns | warns (133/200 in one module) |

If your machine reproduces this, the environment is sound.

### R packages

```r
install.packages(c("ggplot2","dplyr","tidyr","pheatmap","RColorBrewer",
                   "scales","stringr","igraph"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("hgu133plus2.db","AnnotationDbi",
                       "clusterProfiler","org.Hs.eg.db","enrichplot"))
```

---

## Key outputs

| File | Contents |
|---|---|
| `results/FINAL_BIOMARKER_REPORT.md` | **auto-generated report, numbers pre-filled** |
| `results/09_cluster_agreement.csv` | ARI/NMI, clusters vs diagnosis |
| `results/09_cluster_vs_diagnosis.csv` | contingency table |
| `results/09_covariate_association.csv` | **batch check** — what explains each PC |
| `results/05_consensus_pac.csv` | cluster stability |
| `results/08_network_edges_{global,within}.csv` | correlation network |
| `results/enrichment_top.csv` | GO/KEGG, both tests |
| `results/12_candidate_biomarkers.csv` | statistical candidates |
| `results/consensus_biomarkers.csv` | **the deliverable shortlist** |
| `results/14_validation_tiers.csv` | validation strength ladder |
| `figures/fig10_marker_heatmap.png` | **the deliverable heatmap** |

Step 25 writes the report with 9 `TODO` markers — the biology, interpretation and
dataset-specific limitations. Everything countable is already filled in.

---

## Interpreting three steps that mislead people

**Step 14, the correlation network.** In a cohort spanning many subtypes,
gene-gene correlation is dominated by subtype identity: any two genes high in the
same subtype correlate whether or not they are functionally related. You get a few
enormous modules that merely re-describe your clusters. The script computes a
second, within-cluster network that removes the between-subtype signal, and warns
you when one module swallows over half the genes. Report both; interpret the
within-cluster one. WGCNA is the proper tool if you have time.

**Step 24, the ROC curves.** Computed on the held-out split only, using the same
split as step 21, so `--test-size` must match between `s13` and `s14`. A consensus
gene with an AUC near 0.5 that the model still relies on is not a contradiction —
it may be informative only in combination. Say so rather than dropping it.

**Step 25, validation.** Three tiers, weakest to strongest: held-out single genes,
held-out gene panel, independent cohort. Only the third is validation in the
strict sense. `14_validation_tiers.csv` leaves tier 3 blank until you run
`s16_external_validation.py` against a second GEO series.

---

## Known limitations to state in the report

- **Steps 17–20 remain partly circular** if you chose `--group-source cluster`.
  With the default `diagnosis`, they are not.
- **Batch association is measured but not corrected.** No ComBat, no covariate
  term. Step 15 tells you the size of the problem.
- **One-vs-rest testing** makes the reference group heterogeneous, inflating
  significance for genes absent from one large competing subtype.
- **Probe collapsing by max-mean** is one rule among several; a different rule
  gives a slightly different gene set.
- **Bulk tissue** — differences may reflect cell-type proportions rather than
  regulation within a cell type.
- **Without step 25 tier 3, nothing is externally validated.**

## Useful variations

```bash
python3 s04_filter_normalize.py --top-var 2000 --drop-outliers
python3 s05_hierarchical_clustering.py --k 18        # force k = number of diagnoses
python3 s03_quality_assessment.py --skip-qnorm       # ablation for the report
python3 s08_gene_network.py --threshold 0.85         # sparser network
python3 s10_differential_expression.py --fdr 0.01 --lfc 1.5
python3 s12_candidate_biomarkers.py --min-specificity 2.0
python3 s09_metadata_groups.py --group-source cluster  # circular variant, warns
```

Running step 5 twice — once with automatic k, once with `--k 18` — gives the
strongest discussion section. Silhouette almost always prefers fewer clusters than
there are real disease entities, because several subtypes share expression
programmes. That gap is a finding.

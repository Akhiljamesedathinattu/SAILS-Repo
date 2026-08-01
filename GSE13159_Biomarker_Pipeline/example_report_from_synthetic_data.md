# Candidate Biomarkers in Acute Leukaemia

**Dataset** GSE13159 (MILE study) · **Author** Test Run · **Generated** 2026-07-27

> Auto-generated from pipeline outputs. Numbers are final; every **TODO** marks prose you must write yourself.


## 1. Dataset

- Platform: GPL570
- Samples: 180 · Probes: 3000
- Value scale as downloaded: linear (MAS5-like)
- After filtering: 1500 probes retained from 3000


**TODO** one paragraph on why this cohort suits the question.


## 2. Preprocessing and quality control

- Missing values: 0.0 (0.0% of the matrix)
- Mean inter-array correlation after normalisation: 0.008
- Samples flagged as QC outliers: 0 (flagged, not removed)
- log2 transform and quantile normalisation applied within quality assessment, because the data ship on a linear scale and the cohort spans multiple laboratories.


## 3. Clustering

- k = 5 chosen by maximum mean silhouette; lowest PAC would suggest k = 5
- Silhouette 0.026 · mean consensus 0.986


Agreement with curated diagnoses (clusters never saw these labels):

| method | k | adjusted_rand_index | normalized_mutual_info |
|---|---|---|---|
| hclust_cluster | 5 | 0.986 | 0.984 |
| kmeans_cluster | 5 | 1.000 | 1.000 |


Cluster composition:

| hclust_cluster | n | dominant_diagnosis | purity |
|---|---|---|---|
| 1 | 36 | CML | 1.000 |
| 2 | 37 | ALL with t(12;21) | 0.973 |
| 3 | 36 | AML with t(15;17) | 1.000 |
| 4 | 36 | non-leukemia bone marrow | 1.000 |
| 5 | 35 | CLL | 1.000 |


**TODO** which clusters are pure, which fuse two diagnoses, which diagnosis splits across clusters, and whether those merges are biologically sensible.


### Batch check

| covariate | pc | n_levels | eta_squared | p_value |
|---|---|---|---|---|
| leukemia_class | 1 | 5 | 0.907 | 3.95e-34 |
| contact_institute | 1 | 3 | 0.000 | 0.972 |


**TODO** if a technical field rivals diagnosis on PC1, say so explicitly and discuss the confound.


## 4. Gene correlation network

- 200 genes, |r| >= 0.6
- Global edges: 0 · within-cluster edges: 0
- 8 modules, largest holding 133 genes


**TODO** state which network you interpret and why. In a multi-subtype cohort the global network is dominated by subtype identity; the within-cluster network is the more meaningful one.


## 5. Differential expression

| group | n_samples | n_up | n_down |
|---|---|---|---|
| ALL with t(12;21) | 36 | 40 | 22 |
| AML with t(15;17) | 36 | 41 | 28 |
| CLL | 36 | 40 | 25 |
| CML | 36 | 42 | 26 |
| non-leukemia bone marrow | 36 | 41 | 30 |


**TODO** comment on the relation between group size and gene counts; small groups have less power but a homogeneous phenotype can still give very small p-values.


## 6. GO and KEGG enrichment

| library | group | n_query_genes | n_ora_significant | n_rank_significant |
|---|---|---|---|---|
| SYNTHETIC_Pathways | ALL with t(12;21) | 40 | 1 | 5 |
| SYNTHETIC_Pathways | AML with t(15;17) | 41 | 2 | 5 |
| SYNTHETIC_Pathways | CLL | 40 | 1 | 5 |
| SYNTHETIC_Pathways | CML | 42 | 2 | 5 |
| SYNTHETIC_Pathways | non-leukemia bone marrow | 41 | 1 | 5 |


Top terms:

| library | group | label | test | effect | fdr |
|---|---|---|---|---|---|
| SYNTHETIC_Pathways | ALL with t(12;21) | true markers ALL_with_t12;21 | ora | 75.000 | 8.71e-92 |
| SYNTHETIC_Pathways | AML with t(15;17) | true markers AML_with_t15;17 | ora | 73.171 | 1.07e-89 |
| SYNTHETIC_Pathways | AML with t(15;17) | decoy set 13 | ora | 7.317 | 0.011 |
| SYNTHETIC_Pathways | CLL | true markers CLL | ora | 75.000 | 8.71e-92 |
| SYNTHETIC_Pathways | CML | true markers CML | ora | 71.429 | 1.50e-88 |
| SYNTHETIC_Pathways | CML | decoy set 18 | ora | 7.143 | 0.008 |
| SYNTHETIC_Pathways | non-leukemia bone marrow | true markers non-leukemia_bone_marrow | ora | 73.171 | 3.57e-90 |
| SYNTHETIC_Pathways | ALL with t(12;21) | true markers ALL_with_t12;21 | rank | 1.000 | 3.61e-26 |
| SYNTHETIC_Pathways | AML with t(15;17) | true markers AML_with_t15;17 | rank | 1.000 | 3.61e-26 |
| SYNTHETIC_Pathways | CLL | true markers CLL | rank | 1.000 | 3.61e-26 |
| SYNTHETIC_Pathways | CML | true markers CML | rank | 1.000 | 3.61e-26 |
| SYNTHETIC_Pathways | non-leukemia bone marrow | true markers non-leukemia_bone_marrow | rank | 1.000 | 3.61e-26 |


**TODO** lead with terms significant under BOTH tests. State that the background was the filtered gene universe, not the genome.


## 7. Candidate biomarkers

| rank | gene | group | log2FC | fdr | specificity | candidate_score |
|---|---|---|---|---|---|---|
| 1 | GENE00175 | non-leukemia bone marrow | 3.500 | 4.78e-14 | 35.004 | 0.803 |
| 2 | GENE00193 | non-leukemia bone marrow | 3.384 | 1.03e-13 | 11.229 | 0.767 |
| 3 | GENE00085 | CLL | 3.294 | 7.67e-15 | 9.381 | 0.764 |
| 4 | GENE00194 | non-leukemia bone marrow | 3.512 | 6.59e-13 | 12.627 | 0.760 |
| 5 | GENE00042 | AML with t(15;17) | 3.618 | 3.48e-17 | 6.059 | 0.756 |
| 6 | GENE00082 | CLL | 3.146 | 1.26e-14 | 10.935 | 0.756 |
| 7 | GENE00122 | CML | 3.400 | 3.44e-13 | 10.430 | 0.753 |
| 8 | GENE00171 | non-leukemia bone marrow | 3.486 | 2.12e-15 | 7.979 | 0.749 |
| 9 | GENE00155 | CML | 3.465 | 8.94e-14 | 9.137 | 0.742 |
| 10 | GENE00078 | AML with t(15;17) | 3.396 | 9.90e-13 | 33.957 | 0.738 |
| 11 | GENE00149 | CML | 3.476 | 3.44e-13 | 9.395 | 0.735 |
| 12 | GENE00072 | AML with t(15;17) | 3.732 | 1.91e-14 | 7.432 | 0.732 |
| 13 | GENE00052 | AML with t(15;17) | 3.250 | 3.84e-13 | 21.725 | 0.727 |
| 14 | GENE00018 | ALL with t(12;21) | 3.182 | 2.82e-13 | 12.717 | 0.723 |
| 15 | GENE00113 | CLL | 3.230 | 2.42e-13 | 10.313 | 0.723 |
| 16 | GENE00029 | ALL with t(12;21) | 3.386 | 1.78e-12 | 15.460 | 0.720 |
| 17 | GENE00005 | ALL with t(12;21) | 3.251 | 6.46e-13 | 32.513 | 0.718 |
| 18 | GENE00091 | CLL | 3.069 | 6.22e-14 | 11.506 | 0.717 |
| 19 | GENE00174 | non-leukemia bone marrow | 3.179 | 2.04e-13 | 14.808 | 0.717 |
| 20 | GENE00101 | CLL | 3.064 | 1.55e-13 | 10.357 | 0.699 |


## 8. Classification

| model | cv_balanced_accuracy_mean | test_balanced_accuracy | test_macro_auc_ovr | label_source |
|---|---|---|---|---|
| logistic_regression | 1.000 | 1.000 | 1.000 | diagnosis |
| random_forest | 1.000 | 1.000 | 1.000 | diagnosis |


Labels were the curated clinical diagnoses, not the clusters, so these figures are not circular. Feature selection sat inside the cross-validation pipeline and never saw the test split.


## 9. Consensus biomarkers

40 genes carry both statistical and predictive support.

| rank | gene | target_group | log2FC | fdr | importance_rank | univariate_auc_test |
|---|---|---|---|---|---|---|
| 1 | GENE00082 | CLL | 3.146 | 1.26e-14 | 7 | 0.870 |
| 2 | GENE00072 | AML with t(15;17) | 3.732 | 1.91e-14 | 9 | 0.995 |
| 3 | GENE00029 | ALL with t(12;21) | 3.386 | 1.78e-12 | 8 | 0.856 |
| 4 | GENE00193 | non-leukemia bone marrow | 3.384 | 1.03e-13 | 32 | 0.911 |
| 5 | GENE00171 | non-leukemia bone marrow | 3.486 | 2.12e-15 | 25 | 0.973 |
| 6 | GENE00052 | AML with t(15;17) | 3.250 | 3.84e-13 | 21 | 0.916 |
| 7 | GENE00149 | CML | 3.476 | 3.44e-13 | 24 | 0.888 |
| 8 | GENE00042 | AML with t(15;17) | 3.618 | 3.48e-17 | 40 | 0.970 |
| 9 | GENE00125 | CML | 3.361 | 6.25e-12 | 20 | 0.903 |
| 10 | GENE00175 | non-leukemia bone marrow | 3.500 | 4.78e-14 | 49 | 0.962 |
| 11 | GENE00078 | AML with t(15;17) | 3.396 | 9.90e-13 | 39 | 0.959 |
| 12 | GENE00046 | AML with t(15;17) | 3.692 | 5.07e-16 | 2 | 0.959 |
| 13 | GENE00102 | CLL | 3.649 | 1.34e-13 | 19 | 0.979 |
| 14 | GENE00035 | ALL with t(12;21) | 3.223 | 4.13e-12 | 22 | 0.934 |
| 15 | GENE00039 | ALL with t(12;21) | 3.725 | 3.16e-15 | 3 | 0.973 |
| 16 | GENE00122 | CML | 3.400 | 3.44e-13 | 61 | 0.922 |
| 17 | GENE00083 | CLL | 3.268 | 1.05e-15 | 37 | 0.960 |
| 18 | GENE00096 | CLL | 3.364 | 1.05e-15 | 16 | 0.925 |
| 19 | GENE00091 | CLL | 3.069 | 6.22e-14 | 51 | 0.892 |
| 20 | GENE00028 | ALL with t(12;21) | 2.987 | 1.95e-12 | 15 | 0.793 |


Median single-gene AUC: 0.932


**TODO** for each of the top genes, one or two sentences of literature context. Genes already known to mark that subtype are your positive control that the pipeline works; unexpected genes are hypotheses, not discoveries.


## 10. Validation

- 25-gene panel on held-out samples: macro AUC 1.000, balanced accuracy 1.000

| tier | description | metric | value |
|---|---|---|---|
| 1 | single genes, held-out split | median univariate AUC | 0.933 |
| 2 | 25-gene panel, held-out split | macro AUC | 1.000 |
| 2 | all-gene model, held-out split | macro AUC | 1.000 |
| 3 | independent cohort | macro AUC |  |


**No independent cohort was tested.** All figures above are internal to GSE13159. This is a limitation, not a failure, but it must be stated plainly: nothing here is externally validated.


## 11. Limitations

- One-vs-rest testing makes the reference group heterogeneous, inflating significance for genes absent from one large competing subtype.
- Multiple probes per gene were collapsed by max-mean; a different rule would give a slightly different gene set.
- Bulk tissue: differences may reflect changing cell-type proportions rather than regulation within a cell type.
- No batch covariate was modelled; centre association is measured but not corrected.
- **TODO** add any dataset-specific limitations you encountered.


## Figures

_No figures found. Run the R scripts._

## Reproducibility

```bash
export SAILS_BASE=/tmp/bm_test
bash run_all.sh
```

Random seed 42 throughout.

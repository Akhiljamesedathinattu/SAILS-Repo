#!/usr/bin/env bash
# GSE13159 biomarker pipeline — runs all 25 steps in the specified order.
#
#   export SAILS_BASE=/sails/SAILS-Repo/Gene_Expression_Clustering
#   bash run_all.sh
#
# Every step caches its output, so later steps can be re-run without re-parsing
# the download. Safe to interrupt and resume.
set -euo pipefail

export SAILS_BASE="${SAILS_BASE:-/sails/SAILS-Repo/Gene_Expression_Clustering}"
cd "$(dirname "$0")"
mkdir -p "$SAILS_BASE"/{raw,work,results,figures,models}
S() { echo; echo "=== $* ==="; }

echo "project root: $SAILS_BASE"

S "steps 1-2   data understanding"
python3 s01_data_understanding.py

S "step 3      expression matrix loading"
python3 s02_expression_matrix.py

S "steps 4-5   quality assessment + missing values"
python3 s03_quality_assessment.py

S "steps 6-7   variance filtering + z-score"
python3 s04_filter_normalize.py --top-var 5000

S "steps 8-10  hierarchical clustering, dendrogram, cluster assignment"
python3 s05_hierarchical_clustering.py --kmax 20 --consensus-kmax 12

S "step 11     PCA visualisation"
python3 s06_pca.py

S "steps 12-13 probe-to-gene mapping, representative genes"
python3 s07_probe_gene_mapping.py

S "step 14     gene correlation network"
python3 s08_gene_network.py --n-genes 300 --threshold 0.7

S "steps 15-16 clinical metadata integration, patient group creation"
python3 s09_metadata_groups.py --label-field leukemia_class --group-source diagnosis

S "step 17     differential expression"
python3 s10_differential_expression.py

S "steps 18-19 GO + KEGG enrichment"
python3 s11_go_kegg.py || echo "enrichment skipped (no gene sets — see fetch_genesets.py)"

S "step 20     candidate biomarker identification"
python3 s12_candidate_biomarkers.py

S "steps 21-22 machine learning + SHAP"
python3 s13_machine_learning.py --shap

S "steps 23-25 consensus biomarkers, ROC, validation"
python3 s14_consensus_roc.py

S "figures"
if command -v Rscript >/dev/null 2>&1; then
  Rscript figures_core.R
  Rscript figures_ml.R
  Rscript figures_network.R    || echo "network figures skipped"
  Rscript figures_enrichment.R || echo "enrichment figures skipped"
else
  echo "Rscript not found — install R for figures"
fi

S "step 25     final biomarker report"
python3 s15_final_report.py

echo
echo "report : $SAILS_BASE/results/FINAL_BIOMARKER_REPORT.md"
echo "tables : $SAILS_BASE/results"
echo "figures: $SAILS_BASE/figures"
echo
echo "For true external validation (tier 3):"
echo "  python3 s16_external_validation.py --input raw/<other>_series_matrix.txt.gz"

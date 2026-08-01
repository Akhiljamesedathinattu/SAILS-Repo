#!/usr/bin/env Rscript
# Optional but STRONGLY recommended: build raw/GPL570_annot.csv so the pipeline
# works at gene level instead of probe level. Without it, enrichment cannot run
# and the shortlist contains unreadable probe IDs like 205548_s_at.
#
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("hgu133plus2.db", "AnnotationDbi"))
#   Rscript optional_annotate_probes.R
#
# Run this BEFORE 03_collapse_filter_pca.py.

suppressPackageStartupMessages({
  library(hgu133plus2.db); library(AnnotationDbi)
})

BASE <- Sys.getenv("SAILS_BASE", "/home/sails/SAILS-Repo/Gene_Expression_Clustering")
dir.create(file.path(BASE, "raw"), showWarnings = FALSE, recursive = TRUE)
out <- file.path(BASE, "raw", "GPL570_annot.csv")

probes <- keys(hgu133plus2.db, keytype = "PROBEID")
ann <- AnnotationDbi::select(hgu133plus2.db, keys = probes,
                             columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                             keytype = "PROBEID")
ann <- ann[!duplicated(ann$PROBEID), ]
ann <- ann[!is.na(ann$SYMBOL), ]
write.csv(data.frame(probe = ann$PROBEID, symbol = ann$SYMBOL,
                     gene_name = ann$GENENAME, entrez = ann$ENTREZID),
          out, row.names = FALSE)
message("wrote ", out, " (", nrow(ann), " probes with a gene symbol)")

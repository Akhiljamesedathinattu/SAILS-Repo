#!/usr/bin/env Rscript

# =============================================================
# MAKE THE PROBE-TO-GENE FILE
#
# Run this ONCE, BEFORE simple_pipeline.py.
#
# WHY WE NEED IT
# The chip measures things called probes, with codes like
# "205548_s_at". Nobody can read that, and pathway databases
# do not know those codes either. They only know gene names
# like "CD3D".
#
# This script downloads the official list that says which gene
# each probe measures, and saves it as a simple CSV file.
#
# FIRST TIME SETUP (only needed once):
#   install.packages("BiocManager")
#   BiocManager::install(c("hgu133plus2.db", "AnnotationDbi"))
#
# THEN RUN:
#   Rscript simple_make_annotation.R
# =============================================================

suppressPackageStartupMessages(library(hgu133plus2.db))
suppressPackageStartupMessages(library(AnnotationDbi))


# Use the FULL path to your project folder, the same one you put in
# simple_pipeline.py.

BASE_FOLDER <- "/home/sails/SAILS-Repo/Gene_Expression_Clustering"
RAW_FOLDER <- file.path(BASE_FOLDER, "raw")


make_annotation <- function() {
  # Make the raw folder if it is not there yet
  dir.create(RAW_FOLDER, showWarnings = FALSE, recursive = TRUE)

  output_file <- file.path(RAW_FOLDER, "GPL570_annot.csv")

  cat("Getting the list of all probes...\n")
  all_probes <- keys(hgu133plus2.db, keytype = "PROBEID")
  cat("Found", length(all_probes), "probes\n")

  cat("Looking up the gene name for each probe...\n")
  lookup <- AnnotationDbi::select(hgu133plus2.db,
                                  keys = all_probes,
                                  columns = "SYMBOL",
                                  keytype = "PROBEID")

  # Some probes appear more than once. Keep only the first time
  # we see each probe, so every probe has exactly one gene.
  lookup <- lookup[!duplicated(lookup$PROBEID), ]

  # Some probes have no gene name at all. Drop those rows.
  lookup <- lookup[!is.na(lookup$SYMBOL), ]

  # Build a simple two-column table with the names Python expects
  result <- data.frame(probe = lookup$PROBEID,
                       symbol = lookup$SYMBOL)

  write.csv(result, output_file, row.names = FALSE)

  cat("Wrote", output_file, "with", nrow(result), "probes\n")
}


mymain <- function() {
  make_annotation()
}


mymain()

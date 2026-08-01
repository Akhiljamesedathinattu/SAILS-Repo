#!/usr/bin/env Rscript
# =============================================================================
# OPTIONAL — GO and KEGG via clusterProfiler
#
# The Python route (06_enrichment.py) is self-contained and works offline from
# GMT files. This is the Bioconductor alternative, and it gives you things the
# GMT route cannot:
#
#   - real GO ontology structure, so redundant parent/child terms can be pruned
#     with simplify()
#   - live KEGG pathway definitions with pathway IDs, not a frozen snapshot
#   - GSEA on the full ranked list, not just the significant subset
#   - category-network and enrichment-map plots
#
# Run either route, or both and compare — agreement across two independent
# implementations is a genuinely strong result to report.
#
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("clusterProfiler","org.Hs.eg.db","enrichplot","DOSE"))
#   Rscript optional_enrichment_clusterprofiler.R
#
# KEGG queries hit the KEGG REST API, so this script needs internet access.
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(ggplot2); library(dplyr)
})

BASE <- Sys.getenv("SAILS_BASE", "/home/sails/SAILS-Repo/Gene_Expression_Clustering")
RES   <- file.path(BASE, "results")
FIG   <- file.path(BASE, "figures")
GROUP <- Sys.getenv("SAILS_GROUP", "kmeans_cluster")
FDR_IN <- 0.05
LFC_IN <- 1.0
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

de_file <- file.path(RES, paste0("de_results_", GROUP, ".csv"))
if (!file.exists(de_file)) stop("run 05_differential_expression.py first")
de <- read.csv(de_file)

if (!any(grepl("^[A-Z][A-Z0-9-]*$", head(unique(de$feature), 50)))) {
  warning("features do not look like gene symbols — run optional_annotate_probes.R ",
          "and repeat step 03, otherwise nothing will map to Entrez IDs")
}

# ---- symbol -> Entrez, once, for the whole tested universe -------------------
uni_map <- suppressWarnings(
  bitr(unique(as.character(de$feature)), fromType = "SYMBOL",
       toType = "ENTREZID", OrgDb = org.Hs.eg.db))
message(nrow(uni_map), " of ", length(unique(de$feature)),
        " features mapped to Entrez IDs")
universe <- unique(uni_map$ENTREZID)

sym2ent <- setNames(uni_map$ENTREZID, uni_map$SYMBOL)

groups <- sort(unique(de$group))
gene_lists <- lapply(groups, function(g) {
  s <- de %>% filter(group == g, fdr < FDR_IN, log2FC >= LFC_IN) %>% pull(feature)
  unique(na.omit(sym2ent[as.character(s)]))
})
names(gene_lists) <- paste0("C", groups)
gene_lists <- gene_lists[lengths(gene_lists) >= 10]
message("groups with >= 10 mapped DE genes: ", length(gene_lists))
if (length(gene_lists) == 0) stop("no group has enough mapped DE genes")

save_tab <- function(obj, name) {
  if (is.null(obj)) return(invisible(NULL))
  df <- as.data.frame(obj)
  if (nrow(df) == 0) { message("  ", name, ": no significant terms"); return(invisible(NULL)) }
  write.csv(df, file.path(RES, paste0("cp_", name, ".csv")), row.names = FALSE)
  message("  ", name, ": ", nrow(df), " terms -> results/cp_", name, ".csv")
}

# ---- 1. GO over-representation, all three ontologies, compared across groups --
for (ont in c("BP", "MF", "CC")) {
  message("GO:", ont, " compareCluster ...")
  cc <- tryCatch(
    compareCluster(gene_lists, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                   ont = ont, universe = universe, keyType = "ENTREZID",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                   readable = TRUE),
    error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
  if (is.null(cc) || nrow(as.data.frame(cc)) == 0) next

  # collapse redundant parent/child terms — the main advantage over GMT files
  cc_s <- tryCatch(simplify(cc, cutoff = 0.7, by = "p.adjust", select_fun = min),
                   error = function(e) cc)
  save_tab(cc_s, paste0("go_", tolower(ont)))

  p <- dotplot(cc_s, showCategory = 6, label_format = 50) +
    ggtitle(paste0("GO ", ont, " over-representation by cluster")) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(FIG, paste0("fig_cp_go_", tolower(ont), "_dotplot.png")),
         p, width = 12, height = 10, dpi = 300, limitsize = FALSE)
  message("  wrote fig_cp_go_", tolower(ont), "_dotplot.png")
}

# ---- 2. KEGG pathways --------------------------------------------------------
message("KEGG compareCluster (needs internet) ...")
kk <- tryCatch(
  compareCluster(gene_lists, fun = "enrichKEGG", organism = "hsa",
                 universe = universe, pAdjustMethod = "BH", pvalueCutoff = 0.05),
  error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
  kk <- setReadable(kk, org.Hs.eg.db, keyType = "ENTREZID")
  save_tab(kk, "kegg")
  p <- dotplot(kk, showCategory = 8, label_format = 50) +
    ggtitle("KEGG pathway over-representation by cluster") +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(FIG, "fig_cp_kegg_dotplot.png"), p,
         width = 12, height = 9, dpi = 300, limitsize = FALSE)
  message("  wrote fig_cp_kegg_dotplot.png")
}

# ---- 3. GSEA on the full ranked list, per group ------------------------------
for (g in groups) {
  sub <- de %>% filter(group == g) %>%
    mutate(entrez = sym2ent[as.character(feature)]) %>%
    filter(!is.na(entrez)) %>%
    group_by(entrez) %>% slice_max(abs(t_stat), n = 1, with_ties = FALSE) %>% ungroup()
  gl <- setNames(sub$t_stat, sub$entrez)
  gl <- sort(gl, decreasing = TRUE)
  if (length(gl) < 500) next

  message("GSEA cluster ", g, " (", length(gl), " genes) ...")
  gs <- tryCatch(
    gseGO(gl, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID",
          pAdjustMethod = "BH", pvalueCutoff = 0.05, verbose = FALSE),
    error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
  if (is.null(gs) || nrow(as.data.frame(gs)) == 0) next
  save_tab(gs, paste0("gsea_go_bp_cluster", g))

  p <- ridgeplot(gs, showCategory = 15) +
    ggtitle(paste0("GSEA GO:BP — cluster ", g)) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(FIG, paste0("fig_cp_gsea_cluster", g, ".png")),
         p, width = 11, height = 9, dpi = 300, limitsize = FALSE)
}

message("\nclusterProfiler outputs in ", RES, " (cp_*.csv) and ", FIG)
message("Compare cp_kegg.csv against enrichment_ora_KEGG_2021_Human.csv — ",
        "terms found by both routes are the ones worth reporting.")

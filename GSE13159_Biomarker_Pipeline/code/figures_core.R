#!/usr/bin/env Rscript
# =============================================================================
# Core figures (Figures 1-10) — quality control through differential expression
#
# Reads only CSVs produced by the Python steps. No computation on expression
# data happens here.
#
#   install.packages(c("ggplot2","dplyr","tidyr","pheatmap","RColorBrewer","scales"))
#   Rscript figures_core.R
#
# FOUR THINGS THIS FILE HAS TO GET RIGHT, all learned the hard way:
#
# 1. PYTHON BOOLEANS ARE NOT R LOGICALS. pandas writes True/False; read.csv
#    reads those as the character strings "True"/"False", so sum() on them
#    fails with "invalid 'type' (character)". as_logical() below normalises
#    them. Any Python-written boolean column needs this.
#
# 2. THE DATA IS NOT ON A LOG2 SCALE. The GSE13159 series matrix from GEO is
#    min-max scaled to [0,1] per sample. Axis labels here therefore say
#    "normalised expression", and the volcano's x-axis is a difference of
#    normalised means, not a log2 fold change. Hardcoding the usual +/-1
#    fold-change lines would put them off-scale and colour every gene "ns".
#    The cutoff is read from step 10, or derived as a percentile if absent.
#
# 3. SILENT SUCCESS IS WORSE THAN FAILURE. dir.create(showWarnings=FALSE) plus
#    if(!is.null(x)) guards on every plot means a wrong BASE produces a cheerful
#    "figures written" message and zero files. The guards below fail loudly and
#    the summary counts what was actually written.
#
# 4. ASCII ONLY IN pheatmap TITLES. pheatmap renders `main` through grid, whose
#    default device font on this system cannot encode an em dash, producing
#    "conversion failure ... mbcsToSbcs" warnings and substituting dots. ggplot
#    titles are fine; pheatmap titles use a plain hyphen.
#
# CAPTION STYLE: titles are declarative ("Figure 1A. <what it shows>") and the
# subtitle states how to read the panel and what it does or does not establish.
# Written for a thesis figure list, so each caption stands alone without the
# surrounding text.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(pheatmap); library(RColorBrewer); library(scales)
})

BASE <- Sys.getenv("SAILS_BASE", "/home/sails/SAILS-Repo/Gene_Expression_Clustering")
BASE <- normalizePath(path.expand(BASE), mustWork = FALSE)
RES  <- file.path(BASE, "results")
FIG  <- file.path(BASE, "figures")

if (!dir.exists(RES))
  stop("results directory not found: ", RES,
       "\nCheck SAILS_BASE, or run the Python steps first.", call. = FALSE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(FIG))
  stop("could not create figures directory: ", FIG, call. = FALSE)
message("BASE = ", BASE)

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey92", colour = NA),
                  plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9, colour = "grey25")))

.missing <- character(0)
.written <- 0L

rd <- function(f, ...) {
  p <- file.path(RES, f)
  if (!file.exists(p)) { .missing <<- c(.missing, f); return(NULL) }
  read.csv(p, check.names = FALSE, ...)
}
sv <- function(p, f, w = 9, h = 6) {
  ggsave(file.path(FIG, f), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  .written <<- .written + 1L
  message("wrote ", f)
}
noted <- function(f) { .written <<- .written + 1L; message("wrote ", f) }

# Wrap long subtitles so they do not run past the plot edge.
wrap_sub <- function(x, width = 110) paste(strwrap(x, width = width), collapse = "\n")

# Python True/False -> R TRUE/FALSE. Leaves real logicals and 0/1 alone.
as_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

# Effect-size cutoff used by step 10. Written by the corrected s10; otherwise
# recovered from the data so the volcano thresholds match the analysis rather
# than a log2 convention that does not apply to this matrix.
effect_cutoff <- function(vd, col) {
  f <- rd("10_effect_size_cutoff.csv")
  if (!is.null(f) && "cutoff" %in% names(f)) {
    return(list(v = as.numeric(f$cutoff[1]), src = "step 10"))
  }
  list(v = as.numeric(quantile(abs(vd[[col]]), 0.99, na.rm = TRUE)),
       src = "99th percentile of observed differences")
}

# ------------------------------------------------------------------ 1. QC ----
qc <- rd("03_qc_sample_summary.csv")
if (!is.null(qc)) {
  qc$idx <- seq_len(nrow(qc))
  p <- ggplot(qc, aes(idx)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = .35) +
    geom_line(aes(y = median), colour = "steelblue4", linewidth = .3) +
    labs(title = "Figure 1A. Per-Sample Expression Distribution After Quantile Normalisation",
         subtitle = wrap_sub(paste(
           "The shaded band spans the interquartile range and the line marks the median for each array.",
           "Note that quantile normalisation forces every sample to a common empirical distribution,",
           "so a flat profile is expected by construction: this panel confirms that normalisation was",
           "applied and does not by itself demonstrate data quality. Array quality is assessed in Figure 1B.")),
         x = "Sample Index", y = "Normalised Expression")
  sv(p, "fig01a_qc_distribution.png")
  
  if ("mean_correlation" %in% names(qc)) {
    cut_lo <- mean(qc$mean_correlation) - 5 * sd(qc$mean_correlation)
    n_out <- if ("outlier" %in% names(qc)) sum(as_logical(qc$outlier), na.rm = TRUE)
    else sum(qc$mean_correlation < cut_lo, na.rm = TRUE)
    p <- ggplot(qc, aes(mean_correlation)) +
      geom_histogram(bins = 60, fill = "steelblue", colour = "white", linewidth = .2) +
      geom_vline(xintercept = cut_lo, linetype = 2, colour = "#C44E52") +
      labs(title = "Figure 1B. Array Quality Assessed by Mean Inter-Array Correlation",
           subtitle = wrap_sub(sprintf(paste(
             "Each array was correlated against all others and the mean taken. The dashed line marks the",
             "outlier threshold at (mean - 5 SD) = %.3f, beyond which %d of %d samples fall. These samples",
             "were flagged and retained rather than removed, since a rare subtype with an atypical profile",
             "is expected to correlate less well with the cohort."), cut_lo, n_out, nrow(qc))),
           x = "Mean Correlation With All Other Samples", y = "Number of Samples")
    sv(p, "fig01b_qc_correlation.png", h = 5)
  }
}

# ---------------------------------------- 2. covariate association (batch) ----
ca <- rd("09_covariate_association.csv")
if (!is.null(ca) && nrow(ca) > 0) {
  p <- ggplot(ca, aes(factor(pc), reorder(covariate, eta_squared), fill = eta_squared)) +
    geom_tile(colour = "white", linewidth = .4) +
    geom_text(aes(label = ifelse(eta_squared > .05, sprintf("%.2f", eta_squared), "")),
              size = 2.8, colour = "grey15") +
    scale_fill_viridis_c(name = expression(eta^2), option = "C", end = .95) +
    labs(title = "Figure 2. Metadata Variables Explaining Each Principal Component",
         subtitle = wrap_sub(paste(
           "Cell values are eta-squared, the proportion of variance in each component explained by the",
           "metadata field in a one-way analysis of variance. A technical field rivalling diagnosis on PC1",
           "would indicate a batch effect. Only the fields shown were tested; date fields that would",
           "detect submission batch were not, so a batch effect of that kind is not excluded.")),
         x = "Principal Component", y = "Metadata Variable")
  sv(p, "fig02_covariate_association.png", w = 10,
     h = max(4.5, 1.5 + .35 * length(unique(ca$covariate))))
}

# --------------------------------------------------------------- 3. PCA ------
pv <- rd("06_pca_variance.csv")
pc <- rd("09_pca_with_diagnosis.csv")
if (!is.null(pv)) {
  cum3 <- sum(pv$variance_explained[1:3]) * 100
  p <- ggplot(head(pv, 30), aes(factor(pc), variance_explained)) +
    geom_col(fill = "#4C72B0") +
    scale_y_continuous(labels = percent) +
    labs(title = "Figure 3. Variance Explained by Successive Principal Components",
         subtitle = wrap_sub(sprintf(paste(
           "First thirty components of a principal component analysis of the z-scored filtered matrix.",
           "The first three components together explain %.1f%% of total variance, with a pronounced",
           "elbow thereafter."), cum3)),
         x = "Principal Component", y = "Variance Explained")
  sv(p, "fig03_pca_scree.png")
}

if (!is.null(pc) && !is.null(pv)) {
  lab <- sprintf("PC%d (%.1f%%)", 1:3, pv$variance_explained[1:3] * 100)
  # step 09 already merges cluster labels into this file. Joining again creates
  # hclust_cluster.x/.y and the plots below silently vanish, so only join when
  # the column is genuinely absent.
  if (!"hclust_cluster" %in% names(pc)) {
    ca_lab <- rd("05_cluster_assignments.csv")
    if (!is.null(ca_lab))
      pc <- left_join(pc, ca_lab[, c("sample", "hclust_cluster")], by = "sample")
  }
  
  if ("hclust_cluster" %in% names(pc)) {
    p <- ggplot(pc, aes(PC1, PC2, colour = factor(hclust_cluster))) +
      geom_point(size = 1.1, alpha = .75) +
      scale_colour_brewer(palette = "Set1", name = "Cluster") +
      labs(title = "Figure 4A. Principal Component Projection Coloured by Unsupervised Cluster",
           subtitle = wrap_sub(paste(
             "Samples projected onto the first two principal components and coloured by Ward cluster",
             "membership at k = 5. Clusters occupy distinct regions of the projection, although the",
             "boundaries between adjacent clusters are gradual rather than sharply separated.")),
           x = lab[1], y = lab[2])
    sv(p, "fig04a_pca_by_cluster.png")
  }
  
  top_cls <- pc %>% count(diagnosis, sort = TRUE) %>% head(12) %>% pull(diagnosis)
  pc$class_grp <- ifelse(pc$diagnosis %in% top_cls, pc$diagnosis, "other")
  p <- ggplot(pc, aes(PC1, PC2, colour = class_grp)) +
    geom_point(size = 1.1, alpha = .75) +
    scale_colour_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(
      length(unique(pc$class_grp))), name = "Curated Diagnosis") +
    labs(title = "Figure 4B. The Same Projection Coloured by Curated Clinical Diagnosis",
         subtitle = wrap_sub(paste(
           "Identical projection to Figure 4A, recoloured by the independently curated diagnostic label.",
           "Diagnoses were never supplied as input to the dimensionality reduction or the clustering,",
           "so the correspondence between this panel and Figure 4A constitutes validation of the",
           "unsupervised structure rather than a fitted result. The twelve largest classes are shown",
           "individually and the remainder pooled as \"other\".")),
         x = lab[1], y = lab[2]) +
    theme(legend.text = element_text(size = 7))
  sv(p, "fig04b_pca_by_diagnosis.png", w = 11)
  
  if ("PC3" %in% names(pc) && "hclust_cluster" %in% names(pc)) {
    p <- ggplot(pc, aes(PC2, PC3, colour = factor(hclust_cluster))) +
      geom_point(size = 1.1, alpha = .75) +
      scale_colour_brewer(palette = "Set1", name = "Cluster") +
      labs(title = "Figure 4C. Second Versus Third Principal Component by Cluster",
           subtitle = wrap_sub(paste(
             "Cluster separation is weaker on these components than on PC1 and PC2, consistent with",
             "the smaller share of variance they capture.")),
           x = lab[2], y = lab[3])
    sv(p, "fig04c_pca_pc2_pc3.png")
  }
}

# ------------------------------------------------------- 4. choosing k -------
scan <- rd("05_k_selection.csv")
pac  <- rd("05_consensus_pac.csv")
summ <- rd("05_clustering_summary.csv")
if (!is.null(scan)) {
  d <- scan %>% pivot_longer(c(silhouette), names_to = "metric", values_to = "value")
  if (!is.null(pac)) d <- bind_rows(d, transmute(pac, k, metric = "pac", value = pac))
  d$metric <- factor(d$metric, levels = c("silhouette", "pac"),
                     labels = c("Mean Silhouette Width (higher is better)",
                                "Proportion of Ambiguous Clustering (lower is better)"))
  best_sil <- scan$k[which.max(scan$silhouette)]
  best_pac <- if (!is.null(pac)) pac$k[which.min(pac$pac)] else NA_integer_
  
  # Mark the k actually used, which may be neither optimum: the pipeline can be
  # told to prefer stability over compactness, and the figure must not imply the
  # silhouette peak was chosen when it was not.
  k_used <- if (!is.null(summ) && "k_final" %in% names(summ)) summ$k_final[1] else NA
  rule   <- if (!is.null(summ) && "k_rule" %in% names(summ)) summ$k_rule[1] else ""
  
  vl <- data.frame(metric = levels(d$metric),
                   xint = c(best_sil, if (!is.null(pac)) best_pac else best_sil))
  
  sub <- paste0(
    "Dashed red lines mark each criterion's own optimum",
    if (!is.na(best_pac)) sprintf(" (silhouette k = %d; PAC k = %d)", best_sil, best_pac)
    else sprintf(" (silhouette k = %d)", best_sil),
    if (!is.na(k_used)) sprintf("; the solid green line marks k = %d, the value used", k_used) else "",
    if (nzchar(rule)) sprintf(", selected by the rule \"%s\"", rule) else "",
    ". Silhouette rewards compact, well-separated clusters and is weakly discriminating in high",
    " dimensions, where all values are low; PAC measures reproducibility under resampling. Where the",
    " two disagree, the disagreement is reported rather than resolved.")
  
  p <- ggplot(d, aes(k, value)) +
    geom_line(colour = "#4C72B0") + geom_point(size = 1.5, colour = "#4C72B0") +
    geom_vline(data = vl, aes(xintercept = xint), linetype = 2, colour = "#C44E52") +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(title = "Figure 5. Selection of the Number of Clusters by Two Independent Criteria",
         subtitle = wrap_sub(sub, 125),
         x = "Number of Clusters (k)", y = NULL)
  if (!is.na(k_used))
    p <- p + geom_vline(xintercept = k_used, linetype = 1,
                        colour = "#55A868", linewidth = .6)
  sv(p, "fig05_k_selection.png", w = 12, h = 4.8)
}

# ------------------------------------------------- 5. consensus matrix -------
cm <- rd("05_consensus_matrix.csv", row.names = 1)
cma <- rd("05_consensus_matrix_anno.csv")
if (!is.null(cm) && !is.null(cma)) {
  m <- as.matrix(cm)
  ann <- data.frame(Cluster = factor(cma$cluster))
  rownames(ann) <- cma$sample
  ann <- ann[colnames(m), , drop = FALSE]
  nc <- nlevels(ann$Cluster)
  pheatmap(m,
           color = colorRampPalette(c("white", "#4C72B0", "#1A2E4A"))(100),
           annotation_col = ann, annotation_row = ann,
           annotation_colors = list(Cluster = setNames(
             colorRampPalette(brewer.pal(9, "Set1"))(nc), levels(ann$Cluster))),
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = FALSE,
           main = paste("Figure 6. Consensus Matrix: Frequency With Which Each Sample Pair",
                        "Co-Clusters Across Resampling"),
           filename = file.path(FIG, "fig06_consensus_matrix.png"),
           width = 10, height = 9)
  noted("fig06_consensus_matrix.png")
}

# ---------------------------------------------------------- 6. dendrogram ----
mg <- rd("dendro_merge.csv")
if (!is.null(mg)) {
  merge_m <- as.matrix(mg); storage.mode(merge_m) <- "integer"
  hc <- list(merge = merge_m,
             height = rd("dendro_height.csv")$height,
             order  = rd("dendro_order.csv")$order,
             labels = as.character(rd("dendro_labels.csv")$label),
             method = "ward.D2", dist.method = "euclidean")
  class(hc) <- "hclust"
  
  cass <- rd("05_cluster_assignments.csv")
  k_final <- if (!is.null(summ) && "k_final" %in% names(summ)) summ$k_final[1]
  else if (!is.null(cass)) max(cass$hclust_cluster)
  else 4
  # The dendrogram is a random subset of samples, so the full-cohort k may
  # exceed the number of separable branches here; clamp to keep rect.hclust
  # from erroring on a tree with fewer leaves than requested rectangles.
  k_final <- max(2, min(k_final, length(hc$labels) - 1))
  
  png(file.path(FIG, "fig07_dendrogram.png"), width = 2400, height = 1300, res = 200)
  op <- par(mar = c(5, 4, 5, 2))
  plot(hc, labels = FALSE, hang = -1, sub = "",
       xlab = paste0("Samples (random subset of ", length(hc$labels), ")"),
       ylab = "Ward Linkage Height",
       main = paste0("Figure 7. Ward Hierarchical Clustering of the Z-Scored Expression Matrix\n",
                     "Coloured rectangles delimit the k = ", k_final,
                     " clusters used throughout; a random subset is shown for legibility"))
  rect.hclust(hc, k = k_final,
              border = colorRampPalette(brewer.pal(9, "Set1"))(k_final))
  par(op)
  dev.off()
  noted("fig07_dendrogram.png")
}

# ------------------------------------------------ 7. cluster vs diagnosis ----
ct <- rd("09_cluster_vs_diagnosis.csv", row.names = 1)
if (!is.null(ct)) {
  ctm <- as.matrix(ct)
  pheatmap(ctm / pmax(rowSums(ctm), 1),
           color = colorRampPalette(brewer.pal(9, "Blues"))(100),
           display_numbers = ctm, number_format = "%.0f", fontsize_number = 7,
           cluster_rows = TRUE, cluster_cols = FALSE,
           main = paste("Figure 8. Cross-Tabulation of Curated Diagnosis Against Unsupervised",
                        "Cluster (row-normalised; counts shown)"),
           filename = file.path(FIG, "fig08_cluster_vs_diagnosis.png"),
           width = 10.5, height = 10)
  noted("fig08_cluster_vs_diagnosis.png")
}

# ------------------------------------------------------------- 8. volcano ----
vd <- rd("10_volcano_data.csv")
if (!is.null(vd)) {
  # prefer the honestly-named column when the corrected step 10 wrote one
  xcol <- if ("mean_diff" %in% names(vd)) "mean_diff" else "log2FC"
  cut <- effect_cutoff(vd, xcol)
  vd$xv <- vd[[xcol]]
  vd <- vd %>% mutate(negl10p = -log10(pmax(p_value, 1e-300)),
                      status = case_when(fdr < .05 & xv >=  cut$v ~ "Up",
                                         fdr < .05 & xv <= -cut$v ~ "Down",
                                         TRUE ~ "Not significant"))
  vd$status <- factor(vd$status, levels = c("Up", "Down", "Not significant"))
  p <- ggplot(vd, aes(xv, negl10p, colour = status)) +
    geom_point(size = .5, alpha = .5) +
    scale_colour_manual(values = c("Up" = "#C44E52", "Down" = "#4C72B0",
                                   "Not significant" = "grey78"), name = NULL) +
    geom_vline(xintercept = c(-cut$v, cut$v), linetype = 2, linewidth = .3) +
    facet_wrap(~ group, labeller = label_wrap_gen(24)) +
    labs(title = "Figure 9. Differential Expression of Each Diagnostic Group Against All Others",
         subtitle = wrap_sub(sprintf(paste(
           "The horizontal axis is the difference of normalised group means, NOT a log2 fold change:",
           "the source matrix is min-max scaled to [0,1] per sample. Dashed lines mark the effect-size",
           "threshold |difference| >= %s, taken as the %s. Points are coloured at FDR < 0.05 and beyond",
           "that threshold. Welch's t-test, Benjamini-Hochberg adjusted within each comparison."),
           signif(cut$v, 4), cut$src)),
         x = "Difference of Normalised Group Means",
         y = expression(-log[10]~italic(p))) +
    theme(strip.text = element_text(size = 7),
          legend.position = "bottom")
  sv(p, "fig09_volcano.png", w = 12, h = 9)
}

# ------------------------------------------------------ 9. marker heatmap ----
mat <- rd("10_heatmap_matrix.csv", row.names = 1)
hca <- rd("10_heatmap_col_anno.csv")
if (!is.null(mat) && !is.null(hca)) {
  m <- as.matrix(mat)
  rownames(hca) <- hca$sample
  col_ann <- hca[colnames(m), c("group", "cluster"), drop = FALSE]
  names(col_ann) <- c("Diagnosis", "Cluster")
  col_ann$Diagnosis <- factor(col_ann$Diagnosis)
  col_ann$Cluster   <- factor(col_ann$Cluster)
  ng <- nlevels(col_ann$Diagnosis); ncl <- nlevels(col_ann$Cluster)
  pheatmap(m,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
           breaks = seq(-3, 3, length.out = 102),
           annotation_col = col_ann,
           annotation_colors = list(
             Diagnosis = setNames(colorRampPalette(brewer.pal(9, "Set1"))(ng),
                                  levels(col_ann$Diagnosis)),
             Cluster = setNames(colorRampPalette(brewer.pal(9, "Paired"))(ncl),
                                levels(col_ann$Cluster))),
           cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D2",
           show_colnames = FALSE, fontsize_row = 7,
           main = paste("Figure 10. Expression of the Top Markers per Diagnostic Group",
                        "(row z-scored within the displayed subset)"),
           filename = file.path(FIG, "fig10_marker_heatmap.png"),
           width = 13, height = 9)
  noted("fig10_marker_heatmap.png")
}

met <- rd("09_cluster_agreement.csv")
if (!is.null(met)) {
  message("\nClustering agreement with curated diagnoses:")
  print(met, row.names = FALSE)
}

# ------------------------------------------------------------------ summary --
# NOTE: braces are required here. At top level R treats `if (x) expr` as a
# complete statement at the end of the line, so a bare `else` on the next line
# is a parse error. Inside { } it is fine, which is why the else-chains earlier
# in this file work unbraced.
if (length(.missing)) {
  message("\nskipped, inputs not found in ", RES, ":\n  ",
          paste(unique(.missing), collapse = "\n  "))
}
if (.written == 0L) {
  message("\nNO FIGURES WRITTEN - every input CSV was missing. ",
          "Run the Python pipeline steps before this script.")
} else {
  message("\n", .written, " core figures written to ", FIG)
}

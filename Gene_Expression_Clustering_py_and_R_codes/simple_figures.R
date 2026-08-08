#!/usr/bin/env Rscript

# =============================================================
# SIMPLE FIGURES SCRIPT
#
# This R script draws the pictures for our analysis.
#
# IMPORTANT RULE: this script does NO maths on the gene data.
# It only reads the answer files that Python already saved, and
# turns them into pictures. Keeping maths and pictures separate
# means the pictures can never disagree with the numbers.
#
# Run it AFTER simple_pipeline.py:
#     Rscript simple_figures.R
# =============================================================

# ---- Libraries ----
# ggplot2 is the standard R library for making graphs.
# We load it quietly so R does not print startup messages.

suppressPackageStartupMessages(library(ggplot2))


# ---- Settings ----
#
# Use the FULL path to your project folder, exactly the same one you
# put in simple_pipeline.py. A short name like "results" means
# "next to wherever I am standing right now", which breaks as soon as
# you run the script from a different folder.

BASE_FOLDER <- "/home/sails/SAILS-Repo/Gene_Expression_Clustering"

RESULT_FOLDER <- file.path(BASE_FOLDER, "results")
FIGURE_FOLDER <- file.path(BASE_FOLDER, "figures")


make_folder <- function() {
  # Create the figures folder if it does not exist yet.
  # recursive = TRUE lets it make parent folders too.
  # showWarnings = FALSE stops R complaining if it already exists.

  dir.create(FIGURE_FOLDER, showWarnings = FALSE, recursive = TRUE)
  cat("Figures folder ready\n")
}


read_result <- function(file_name) {
  # Read one of the CSV files Python saved.
  # If the file is missing we return NULL instead of crashing,
  # so one missing file does not stop all the other pictures.

  full_path <- file.path(RESULT_FOLDER, file_name)

  if (!file.exists(full_path)) {
    cat("Missing file, skipping:", file_name, "\n")
    return(NULL)
  }

  data <- read.csv(full_path)
  return(data)
}


save_picture <- function(picture, file_name, w = 9, h = 6) {
  # Save a ggplot picture as a PNG image.
  # dpi = 300 makes it sharp enough for a report or thesis.

  full_path <- file.path(FIGURE_FOLDER, file_name)
  ggsave(full_path, picture, width = w, height = h, dpi = 300)
  cat("Wrote", file_name, "\n")
}


draw_quality_plot <- function() {
  # Picture 1: how similar is each chip to all the others?
  # A dot far below the crowd is an unusual chip.

  quality <- read_result("step2_quality.csv")
  if (is.null(quality)) {
    return(invisible(NULL))
  }

  # Add a simple counting column so we have something for the x axis
  quality$chip_number <- 1:nrow(quality)

  picture <- ggplot(quality, aes(x = chip_number, y = mean_correlation)) +
    geom_point(size = 0.6, colour = "steelblue") +
    labs(title = "Chip quality check",
         x = "Chip number",
         y = "Average similarity to other chips") +
    theme_bw()

  save_picture(picture, "figure1_quality.png")
}


draw_k_plot <- function() {
  # Picture 2: how good is each possible number of patient groups?
  # The highest point is the best answer.

  scores <- read_result("step4_tightness.csv")
  if (is.null(scores)) {
    return(invisible(NULL))
  }

  picture <- ggplot(scores, aes(x = k, y = silhouette)) +
    geom_line(colour = "grey40") +
    geom_point(size = 2, colour = "darkred") +
    labs(title = "Choosing the number of patient groups",
         x = "Number of groups (k)",
         y = "Silhouette score (higher is better)") +
    theme_bw()

  save_picture(picture, "figure2_choosing_k.png")
}


draw_pca_plot <- function() {
  # Picture 3: every patient as one dot in 2D.
  # Patients with similar genes sit close together.
  # Colours show which cluster we put them in.

  pca <- read_result("step5_pca.csv")
  if (is.null(pca)) {
    return(invisible(NULL))
  }

  # The cluster number is a number, but we want R to treat it as a
  # category (a label), not a quantity. as.factor does that.
  pca$cluster <- as.factor(pca$tree_cluster)

  picture <- ggplot(pca, aes(x = PC1, y = PC2, colour = cluster)) +
    geom_point(size = 1, alpha = 0.7) +
    labs(title = "Patients grouped by their gene patterns",
         x = "PC1", y = "PC2", colour = "Cluster") +
    theme_bw()

  save_picture(picture, "figure3_pca.png")
}


draw_volcano_plot <- function() {
  # Picture 4: a volcano plot.
  #
  # x axis = how big the difference is
  # y axis = how sure we are the difference is real
  #
  # The interesting genes are in the top corners: big difference
  # AND strong evidence.

  genes <- read_result("step8_volcano.csv")
  cutoff_table <- read_result("step8_cutoff.csv")

  if (is.null(genes) || is.null(cutoff_table)) {
    return(invisible(NULL))
  }

  # Read the cut-off Python used, so our lines match the analysis.
  # We do NOT type a number here by hand.
  size_cutoff <- cutoff_table$cutoff[1]

  # There can be 250,000 rows. Drawing them all is slow, so we take
  # one disease group as an example.
  first_group <- genes$group[1]
  one_group <- genes[genes$group == first_group, ]

  # -log10 turns a tiny p-value into a big number, which is easier
  # to see on a graph. p = 0.001 becomes 3.
  one_group$evidence <- -log10(one_group$p_value + 1e-300)

  picture <- ggplot(one_group, aes(x = difference, y = evidence)) +
    geom_point(size = 0.4, alpha = 0.3, colour = "grey30") +
    geom_vline(xintercept = size_cutoff, colour = "red", linetype = "dashed") +
    geom_vline(xintercept = -size_cutoff, colour = "red", linetype = "dashed") +
    labs(title = paste("Genes that differ in:", first_group),
         x = "Difference in normalised expression",
         y = "Evidence, -log10(p)") +
    theme_bw()

  save_picture(picture, "figure4_volcano.png")
}


draw_importance_plot <- function() {
  # Picture 5: the genes the machine learning model relied on most.

  importance <- read_result("step10_gene_importance.csv")
  if (is.null(importance)) {
    return(invisible(NULL))
  }

  # The file now holds BOTH models, so pick out just one of them.
  # Without this we would mix two different kinds of score on one
  # chart, which would be meaningless.
  #
  # We do NOT type a model name here. Step 10 decided which model to
  # use and wrote that decision to a file, so we read it. Typing a
  # name would risk the picture disagreeing with the analysis.
  chosen <- read_result("step10_chosen_model.csv")

  if (is.null(chosen)) {
    cat("No step10_chosen_model.csv - skipping the importance plot\n")
    return(invisible(NULL))
  }

  wanted_model <- as.character(chosen$value[1])
  importance <- importance[importance$model == wanted_model, ]

  if (nrow(importance) == 0) {
    cat("No rows for model", wanted_model, "- skipping\n")
    return(invisible(NULL))
  }

  # Keep only the top 20 rows
  top_genes <- head(importance, 20)

  # reorder makes the bars sort themselves by size instead of
  # appearing in alphabetical order.
  picture <- ggplot(top_genes,
                    aes(x = reorder(gene, importance), y = importance)) +
    geom_col(fill = "darkgreen") +
    coord_flip() +          # flip so the gene names read left to right
    labs(title = paste("Genes the model used most:", wanted_model),
         x = "Gene", y = "Importance (SHAP)") +
    theme_bw()

  save_picture(picture, "figure5_importance.png")
}


draw_stability_plot <- function() {
  # Picture 6: the stability (PAC) score for each number of groups.
  # LOWER is better here, which is the opposite of the tightness plot,
  # so the axis label says so explicitly.

  scores <- read_result("step4_stability.csv")
  if (is.null(scores)) {
    return(invisible(NULL))
  }

  picture <- ggplot(scores, aes(x = k, y = pac)) +
    geom_line(colour = "grey40") +
    geom_point(size = 2, colour = "darkblue") +
    labs(title = "Choosing the number of groups by stability",
         x = "Number of groups (k)",
         y = "PAC score (LOWER is better)") +
    theme_bw()

  save_picture(picture, "figure6_stability.png")
}


draw_biomarker_roc_plot <- function() {
  # Picture 7: how well each biomarker works as a test.
  #
  # The diagonal dashed line is what pure guessing looks like. A
  # curve that hugs the top-left corner is a good test.

  curves <- read_result("step11_roc_curves.csv")
  if (is.null(curves)) {
    return(invisible(NULL))
  }

  # Put the gene name and its score together in the legend
  curves$label <- paste0(curves$gene, " (", round(curves$auc, 2), ")")

  picture <- ggplot(curves, aes(x = false_positive_rate,
                               y = true_positive_rate,
                               colour = label)) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.8) +
    labs(title = "How well each biomarker works as a test",
         x = "False alarms", y = "Correct catches", colour = "Gene (AUC)") +
    theme_bw()

  save_picture(picture, "figure7_biomarker_roc.png", w = 10, h = 7)
}


draw_network_plot <- function() {
  # Picture 8: which genes have the most connections (the hubs).
  #
  # We show the two networks side by side. The "everyone" bars are
  # inflated by disease identity; the "within" bars are the real
  # gene teamwork. Seeing them together makes the difference obvious.

  genes <- read_result("step12_network_genes.csv")
  if (is.null(genes)) {
    return(invisible(NULL))
  }

  top_genes <- head(genes, 20)

  # Stack the two counts into one long table so ggplot can put them
  # side by side. We build it by hand rather than using tidyr, so
  # there is one less library to install.
  long_table <- data.frame(
    gene = c(top_genes$gene, top_genes$gene),
    network = c(rep("across everyone", nrow(top_genes)),
                rep("within clusters", nrow(top_genes))),
    connections = c(top_genes$lines_everyone, top_genes$lines_within))

  picture <- ggplot(long_table,
                    aes(x = reorder(gene, connections),
                        y = connections, fill = network)) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(title = "Genes with the most connections",
         x = "Gene", y = "Number of connections", fill = "Network") +
    theme_bw()

  save_picture(picture, "figure8_network.png")
}


mymain <- function() {
  cat("Drawing the figures\n")

  make_folder()

  draw_quality_plot()
  draw_k_plot()
  draw_stability_plot()
  draw_pca_plot()
  draw_volcano_plot()
  draw_importance_plot()
  draw_biomarker_roc_plot()
  draw_network_plot()

  cat("Finished. Pictures are in the", FIGURE_FOLDER, "folder\n")
}


mymain()

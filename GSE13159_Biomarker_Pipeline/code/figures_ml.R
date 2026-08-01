#!/usr/bin/env Rscript
# =============================================================================
# Machine learning and biomarker figures (Figures 11-19)
#
#   Rscript figures_ml.R
#
# FOUR ISSUES THIS FILE HAS TO HANDLE:
#
# 1. SCHEMA DRIFT. s14 writes `gene` and `target_group`; the plots below were
#    originally written against an older schema using `feature` and
#    `best_group`. Aliases are added on read rather than renaming downstream.
#
# 2. AMBIGUOUS ROC LEGEND. Both models contribute curves to Figure 12, and a
#    legend keyed on class alone produces duplicate entries with different AUCs
#    (the same class scored by each model) that a reader cannot attribute. The
#    legend key therefore includes the model.
#
# 3. BREWER PALETTES SATURATE. "Set2" supplies eight colours. Any variable with
#    more levels silently drops bars, so palettes are interpolated with
#    colorRampPalette to the observed number of levels.
#
# 4. ASCII ONLY IN pheatmap TITLES. grid's default device font here cannot
#    encode an em dash and emits "conversion failure ... mbcsToSbcs".
#
# CAPTION STYLE matches figures_core.R: declarative title, subtitle explaining
# how to read the panel and what it does or does not establish.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(pheatmap); library(RColorBrewer)
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

wrap_sub <- function(x, width = 110) paste(strwrap(x, width = width), collapse = "\n")

# Interpolate a Brewer palette to however many levels are actually present.
pal_n <- function(n, name = "Set2", k = 8) colorRampPalette(brewer.pal(k, name))(n)

# Human-readable model names for titles and legends.
pretty_model <- function(x) recode(as.character(x),
                                   logistic_regression = "Logistic Regression",
                                   random_forest       = "Random Forest",
                                   .default = as.character(x))

# --------------------------------------------------------- 1. model metrics ----
mm <- rd("ml_metrics.csv")
if (!is.null(mm)) {
  n_feat <- if ("k_features" %in% names(mm)) mm$k_features[1] else NA
  n_tr   <- if ("n_train" %in% names(mm)) mm$n_train[1] else NA
  n_te   <- if ("n_test"  %in% names(mm)) mm$n_test[1]  else NA
  n_cls  <- if ("n_classes" %in% names(mm)) mm$n_classes[1] else NA
  selm   <- if ("selection" %in% names(mm)) as.character(mm$selection[1]) else ""
  
  d <- mm %>%
    dplyr::select(model, cv = cv_balanced_accuracy_mean,
                  test = test_balanced_accuracy, auc = test_macro_auc_ovr) %>%
    mutate(model = pretty_model(model)) %>%
    pivot_longer(-model, names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                           cv   = "Cross-validated balanced accuracy",
                           test = "Held-out balanced accuracy",
                           auc  = "Held-out macro AUC (one-vs-rest)"))
  p <- ggplot(d, aes(model, value, fill = metric)) +
    geom_col(position = "dodge", width = .7) +
    geom_text(aes(label = sprintf("%.3f", value)),
              position = position_dodge(width = .7), vjust = -0.4, size = 3) +
    scale_fill_brewer(palette = "Set2") +
    coord_cartesian(ylim = c(0, 1.08)) +
    labs(title = "Figure 11. Classification Performance for Prediction of Curated Diagnosis",
         subtitle = wrap_sub(sprintf(paste(
           "Class labels are the independently curated clinical diagnoses, not cluster assignments, so",
           "these estimates are not circular. Feature selection was performed inside the",
           "cross-validation pipeline and therefore never saw the held-out split%s. Training n = %s,",
           "held-out n = %s, %s classes. Balanced accuracy is the macro-average of per-class recall,",
           "appropriate given the wide range of class sizes; chance is approximately %.3f."),
           if (nzchar(selm)) sprintf(" (%s feature selection, %s features)", selm, n_feat) else "",
           n_tr, n_te, n_cls, if (is.na(n_cls)) NA_real_ else 1 / as.numeric(n_cls))),
         x = NULL, y = NULL, fill = NULL)
  sv(p, "fig11_ml_metrics.png", w = 10.5, h = 6)
  message("model metrics:"); print(mm, row.names = FALSE)
}

# ------------------------------------------------------------- 2. ROC curves ----
roc <- rd("ml_roc_curves.csv")
if (!is.null(roc)) {
  # Key the legend on model AND class. Keying on class alone yields duplicate
  # entries with different AUCs (one per model) that cannot be attributed.
  roc <- roc %>%
    mutate(model_pretty = pretty_model(model),
           lab = sprintf("%s - %s (%.2f)", model_pretty, class, auc))
  ord <- roc %>% distinct(model_pretty, class, auc, lab) %>%
    arrange(model_pretty, desc(auc)) %>% pull(lab)
  roc$lab <- factor(roc$lab, levels = ord)
  
  p <- ggplot(roc, aes(fpr, tpr, colour = lab)) +
    geom_abline(linetype = 2, colour = "grey60", linewidth = .3) +
    geom_line(linewidth = .5) +
    facet_wrap(~ model_pretty) +
    scale_colour_manual(values = pal_n(nlevels(roc$lab), "Paired", 9),
                        name = "Model - Class (AUC)") +
    coord_equal() +
    labs(title = "Figure 12. One-Versus-Rest ROC Curves on the Held-Out Test Split",
         subtitle = wrap_sub(paste(
           "One curve per diagnostic class per model. Legend entries name the model as well as the",
           "class, because the same class is scored by both models and would otherwise appear twice",
           "with different areas and no way to attribute them. The diagonal marks chance.")),
         x = "False Positive Rate", y = "True Positive Rate") +
    theme(legend.text = element_text(size = 6),
          legend.key.height = unit(8, "pt"))
  sv(p, "fig12_ml_roc.png", w = 14, h = 7)
}

# ------------------------------------------------------- 3. confusion matrix ----
cm_panel <- c(logistic_regression = "13A", random_forest = "13B")
for (mdl in c("logistic_regression", "random_forest")) {
  cm <- rd(paste0("ml_confusion_", mdl, ".csv"), row.names = 1)
  if (is.null(cm)) next
  m <- as.matrix(cm)
  pheatmap(m / pmax(rowSums(m), 1),
           color = colorRampPalette(brewer.pal(9, "Blues"))(100),
           display_numbers = m, number_format = "%.0f", fontsize_number = 7,
           cluster_rows = FALSE, cluster_cols = FALSE,
           main = paste0("Figure ", cm_panel[[mdl]], ". Confusion Matrix, ",
                         pretty_model(mdl),
                         " (rows = true class, row-normalised; counts shown)"),
           filename = file.path(FIG, paste0("fig13_confusion_", mdl, ".png")),
           width = 10.5, height = 9)
  noted(paste0("fig13_confusion_", mdl, ".png"))
}

# --------------------------------------------------------- 4. importance ----
imp <- rd("ml_feature_importance.csv")
if (!is.null(imp)) {
  meth <- paste(unique(imp$method), collapse = ", ")
  d <- imp %>%
    mutate(model = pretty_model(model)) %>%
    group_by(model) %>% slice_max(importance, n = 25) %>% ungroup()
  p <- ggplot(d, aes(importance, reorder(feature, importance))) +
    geom_col(fill = "#4C72B0", width = .7) +
    facet_wrap(~ model, scales = "free") +
    labs(title = "Figure 14. Twenty-Five Most Important Genes per Classifier",
         subtitle = wrap_sub(sprintf(paste(
           "Importance measure: %s. Values are averaged across all one-versus-rest classifiers, so a",
           "gene that sharply identifies a single small subtype is diluted here; per-class attribution",
           "is required to recover such genes and is reported separately. Note that the two models",
           "select largely different genes despite comparable accuracy, which is expected when many",
           "correlated genes carry overlapping information."), meth)),
         x = "Mean Absolute Importance", y = NULL) +
    theme(axis.text.y = element_text(size = 7))
  sv(p, "fig14_feature_importance.png", w = 12.5, h = 7.5)
}

# ----------------------------------------------- 5. consensus biomarkers ----
cb <- rd("consensus_biomarkers.csv")
if (!is.null(cb)) {
  if (!"feature" %in% names(cb) && "gene" %in% names(cb))
    cb$feature <- cb$gene
  if (!"best_group" %in% names(cb) && "target_group" %in% names(cb))
    cb$best_group <- cb$target_group
  
  p <- ggplot(cb, aes(de_percentile, importance_percentile)) +
    geom_point(aes(size = univariate_auc_test, colour = consensus_score),
               alpha = .85) +
    scale_colour_viridis_c(name = "Consensus\nscore", option = "C", end = .92) +
    scale_size_continuous(name = "Single-gene\nAUC", range = c(1.5, 6)) +
    geom_text(data = head(cb, 10), aes(label = feature),
              size = 2.7, vjust = -1.1, colour = "grey20") +
    labs(title = "Figure 15. Statistical Versus Predictive Evidence for Each Candidate Biomarker",
         subtitle = wrap_sub(paste(
           "Each point is one gene. The horizontal axis is its percentile rank by differential",
           "expression, the vertical axis its percentile rank by model importance. Genes towards the",
           "upper right are supported by both lines of evidence; the consensus score is the harmonic",
           "mean of the two percentiles and therefore penalises one-sided support. The ten",
           "highest-scoring genes are labelled. Scatter along the axes reflects redundancy: a gene can",
           "be individually informative yet contribute little once a correlated neighbour is in the model.")),
         x = "Differential Expression Percentile",
         y = "Model Importance Percentile")
  sv(p, "fig15_consensus_biomarkers.png", w = 10.5, h = 7.5)
  
  top <- head(cb, 25)
  ngrp <- length(unique(top$best_group))
  p <- ggplot(top, aes(univariate_auc_test, reorder(feature, univariate_auc_test))) +
    geom_col(aes(fill = factor(best_group)), width = .7) +
    geom_vline(xintercept = .5, linetype = 2, colour = "grey40") +
    scale_fill_manual(values = pal_n(ngrp, "Set2"), name = "Marker of") +
    coord_cartesian(xlim = c(0.4, 1)) +
    labs(title = "Figure 16. Single-Gene Discriminative Ability on Held-Out Samples",
         subtitle = wrap_sub(paste(
           "Area under the one-versus-rest ROC curve for each consensus biomarker considered alone,",
           "computed on the held-out test split. The dashed line marks chance. A gene close to chance",
           "in isolation may still contribute in combination with others, so these values bound",
           "single-marker performance rather than the panel's.")),
         x = "Univariate AUC (held-out test split)", y = NULL) +
    theme(axis.text.y = element_text(size = 8))
  sv(p, "fig16_biomarker_auc.png", w = 10.5, h = 7.5)
}

# --------------------------------------------------- 6. biomarker ROC/boxes ----
br <- rd("biomarker_roc_curves.csv")
if (!is.null(br)) {
  # Order legend by AUC and drop unused levels so the key lists only genes that
  # actually appear in a panel.
  br <- br %>% mutate(lab = sprintf("%s (%.2f)", feature, auc))
  ord <- br %>% distinct(feature, auc, lab) %>% arrange(desc(auc)) %>% pull(lab)
  br$lab <- factor(br$lab, levels = ord)
  br$lab <- droplevels(br$lab)
  
  p <- ggplot(br, aes(fpr, tpr, colour = lab)) +
    geom_abline(linetype = 2, colour = "grey60", linewidth = .3) +
    geom_line(linewidth = .5) +
    facet_wrap(~ group, labeller = label_wrap_gen(26)) + coord_equal() +
    scale_colour_manual(values = pal_n(nlevels(br$lab), "Dark2"),
                        name = "Gene (AUC)", drop = TRUE) +
    labs(title = "Figure 17. Per-Gene ROC Curves for the Consensus Biomarkers",
         subtitle = wrap_sub(paste(
           "Each panel is one diagnostic group and contains only the biomarkers assigned to that group,",
           "so any single gene appears in exactly one panel. Curves are computed on the held-out test",
           "split. The diagonal marks chance.")),
         x = "False Positive Rate", y = "True Positive Rate") +
    theme(legend.text = element_text(size = 7),
          strip.text = element_text(size = 8))
  sv(p, "fig17_biomarker_roc.png", w = 12.5, h = 7.5)
}

be <- rd("biomarker_expression.csv")
if (!is.null(be)) {
  # in_group holds the target group name or "rest" - several values, not two.
  # Collapse to a genuine two-level contrast so the colours and the figure's
  # own title ("target group vs rest") match what is plotted.
  be$side <- ifelse(be$in_group == "rest", "All other samples", "Target group")
  be$side <- factor(be$side, levels = c("Target group", "All other samples"))
  p <- ggplot(be, aes(side, expression, fill = side)) +
    geom_boxplot(outlier.size = .4, width = .6, linewidth = .3) +
    facet_wrap(~ feature, scales = "free_y") +
    scale_fill_manual(values = c("Target group" = "#C44E52",
                                 "All other samples" = "grey75"), guide = "none") +
    labs(title = "Figure 18. Expression of Each Consensus Biomarker in Its Target Group Versus All Others",
         subtitle = wrap_sub(paste(
           "Boxes span the interquartile range with the median marked; whiskers extend to 1.5 times the",
           "interquartile range and points beyond are plotted individually. The vertical axis is",
           "normalised expression on the [0,1] scale of the source matrix, not log2 expression. Each",
           "panel uses an independent vertical scale.")),
         x = NULL, y = "Normalised Expression") +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7),
          strip.text = element_text(size = 8))
  sv(p, "fig18_biomarker_boxplots.png", w = 12.5, h = 8.5)
}

# ------------------------------------------------- 7. external validation ----
em <- rd("external_metrics.csv")
if (!is.null(em) && !is.null(mm)) {
  d <- bind_rows(
    mm %>% transmute(cohort = "Internal held-out split", model,
                     balanced_accuracy = test_balanced_accuracy,
                     macro_auc = test_macro_auc_ovr),
    em %>% transmute(cohort = "Independent external cohort", model,
                     balanced_accuracy = external_balanced_accuracy,
                     macro_auc = external_macro_auc_ovr)) %>%
    mutate(model = pretty_model(model),
           cohort = factor(cohort, levels = c("Internal held-out split",
                                              "Independent external cohort"))) %>%
    pivot_longer(c(balanced_accuracy, macro_auc),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                           balanced_accuracy = "Balanced accuracy",
                           macro_auc = "Macro AUC (one-vs-rest)"))
  p <- ggplot(d, aes(cohort, value, fill = metric)) +
    geom_col(position = "dodge", width = .7) +
    geom_text(aes(label = sprintf("%.3f", value)),
              position = position_dodge(width = .7), vjust = -.4, size = 3) +
    facet_wrap(~ model) +
    scale_fill_brewer(palette = "Set2") +
    coord_cartesian(ylim = c(0, 1.08)) +
    labs(title = "Figure 19. Internal Versus External Cohort Performance",
         subtitle = wrap_sub(paste(
           "The internal estimate is a held-out split of the discovery cohort and shares platform,",
           "processing site and patient population with the training data; the external estimate does",
           "not. A performance drop is therefore expected, and its magnitude is the quantity of",
           "interest rather than a shortcoming to be minimised.")),
         x = NULL, y = NULL, fill = NULL)
  sv(p, "fig19_external_validation.png", w = 11.5, h = 6)
}

# ------------------------------------------------------------------ summary --
if (length(.missing)) {
  message("\nskipped, inputs not found in ", RES, ":\n  ",
          paste(unique(.missing), collapse = "\n  "))
}
if (.written == 0L) {
  message("\nNO FIGURES WRITTEN - every input CSV was missing. ",
          "Run the Python pipeline steps before this script.")
} else {
  message("\n", .written, " ML figures written to ", FIG)
}

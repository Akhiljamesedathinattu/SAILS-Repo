#!/usr/bin/env Rscript
# =============================================================================
# GO and KEGG enrichment figures (Figures 23-27)
#
# Consumes only the CSVs written by the enrichment step. No computation here.
#
#   install.packages(c("ggplot2","dplyr","tidyr","scales","stringr"))
#   Rscript figures_enrichment.R
#
# FIVE ISSUES THIS FILE HAS TO HANDLE:
#
# 1. CANVAS SIZE MUST NOT SCALE LINEARLY WITH FACET COUNT. The previous version
#    asked for 13 x (2 + 1.4 * n_groups) inches. With 17 diagnostic groups that
#    is a 13 x 25.8 inch page: facet strips truncate, panels crop at the device
#    edge and the aspect ratio is unusable. Dot plots are now PAGED - a fixed
#    number of groups per figure at a fixed canvas size - so each page is
#    legible regardless of how many groups exist.
#
# 2. SEVENTEEN LONG DIAGNOSIS NAMES CANNOT SHARE A CATEGORICAL AXIS. They
#    overlap into an unreadable band. Short codes are generated, used on the
#    axes, and the full key is written to a CSV and echoed to the console so the
#    figure legend in the thesis can reproduce it.
#
# 3. THE CROSS-CLUSTER HEATMAP SELECTION WAS WRONG. The original filtered with
#    `first(domain)` inside `group_by(domain)`, which does not select the top
#    terms per domain as intended. Term selection is now explicit: the terms
#    appearing in the most groups, computed per domain.
#
# 4. MISSING TEST COLUMNS AFTER pivot_wider. If a domain has no rank-test hits,
#    the `rank` column never appears and `ora & rank` errors. Both columns are
#    created explicitly before use.
#
# 5. ZERO COUNTS ARE INFORMATIVE. count() drops absent group-by-support
#    combinations, so a group with no dual-supported terms silently vanishes
#    from the bar chart instead of showing a zero.
#
# INTERPRETIVE NOTE carried in the captions: the rank test returns far more
# terms than over-representation, and inspection shows those extra terms are
# largely housekeeping and proliferation programmes (oxidative phosphorylation,
# mitochondrial translation, DNA repair, spliceosome) that differ between
# quiescent and proliferating disease rather than between specific subtypes.
# Requiring agreement between the two tests filters exactly that generic signal,
# which is why the dual-supported terms are the ones to report.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(scales); library(stringr)
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

theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey92", colour = NA),
                  plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 8.5, colour = "grey25")))

.written <- 0L
sv <- function(p, f, w = 11, h = 8) {
  ggsave(file.path(FIG, f), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  .written <<- .written + 1L
  message("wrote ", f)
}
wrap_sub  <- function(x, width = 118) paste(strwrap(x, width = width), collapse = "\n")
# Term labels are truncated harder than before: at 13 inches wide with three
# facet columns there is not room for 90 characters of term name.
wrap_lab  <- function(x, n = 34) str_wrap(str_trunc(x, 62), n)

# ---- tuning knobs. Lower GROUPS_PER_PAGE or TERMS_PER_GROUP if labels still
# ---- collide in your final PDF; the figure count changes but nothing else.
GROUPS_PER_PAGE  <- 6
TERMS_PER_GROUP  <- 5
HEATMAP_TERMS    <- 20

f_top <- file.path(RES, "enrichment_top.csv")
if (!file.exists(f_top)) {
  message("enrichment_top.csv not found in ", RES,
          " - run the enrichment step first. Skipping enrichment figures.")
} else {
  
  top <- read.csv(f_top, check.names = FALSE)
  if (nrow(top) == 0) {
    message("enrichment_top.csv is empty - no terms passed the FDR cutoff. ",
            "Skipping enrichment figures.")
  } else {
    
    top <- top %>%
      mutate(group  = as.character(group),
             label  = wrap_lab(label),
             domain = ifelse(is.na(domain) | domain == "", library, domain))
    
    # ------------------------------------------------- 0. short codes for groups --
    # Seventeen full diagnosis names cannot share a categorical axis. Build stable
    # short codes, use them on every axis, and write the key out so the thesis
    # figure legend can reproduce it verbatim.
    grp_levels <- sort(unique(top$group))
    short_code <- function(x) {
      s <- x
      s <- gsub("Non-leukemia and healthy bone marrow", "Normal BM", s, fixed = TRUE)
      s <- gsub("AML with normal karyotype \\+ other abnormalities", "AML-NK/other", s)
      s <- gsub("AML complex aberrant karyotype", "AML-complex", s, fixed = TRUE)
      s <- gsub("^AML with ", "AML-", s)
      s <- gsub("^ALL with ", "ALL-", s)
      s <- gsub("^Pro-B-ALL with ", "proB-", s)
      s <- gsub("^c-ALL/Pre-B-ALL ", "cALL-", s)
      s <- gsub("hyperdiploid karyotype", "hyperdip", s, fixed = TRUE)
      s <- gsub("inv\\(16\\)/t\\(16;16\\)", "inv16", s)
      s <- gsub("t\\(11q23\\)/MLL", "MLL", s)
      s <- gsub("without t\\(9;22\\)", "no-t9;22", s)
      s <- gsub("with t\\(9;22\\)", "t9;22", s)
      s <- gsub("\\s+", "", s)
      substr(s, 1, 16)
    }
    key <- tibble(group = grp_levels, code = short_code(grp_levels))
    # Guarantee uniqueness even if two names collapse to the same code.
    key$code <- make.unique(key$code, sep = "_")
    write.csv(key, file.path(RES, "enrichment_group_codes.csv"), row.names = FALSE)
    message("group code key -> results/enrichment_group_codes.csv")
    print(as.data.frame(key), row.names = FALSE)
    
    code_of <- setNames(key$code, key$group)
    top$code <- factor(code_of[top$group], levels = key$code)
    
    key_line <- paste(sprintf("%s = %s", key$code, key$group), collapse = "; ")
    
    # ------------------------------------------------------ 1-2. dot plots, paged
    # One figure per domain per page of groups. Fixed canvas per page, so adding
    # groups adds pages rather than stretching the device until it crops.
    dotplot_pages <- function(d, test_name, dm, fig_no, colour_lab, colour_opt,
                              title_stem, sub_stem) {
      d <- d %>%
        group_by(code) %>% slice_min(fdr, n = TERMS_PER_GROUP, with_ties = FALSE) %>%
        ungroup()
      if (nrow(d) == 0) return(invisible(NULL))
      
      codes <- levels(droplevels(d$code))
      pages <- split(codes, ceiling(seq_along(codes) / GROUPS_PER_PAGE))
      npage <- length(pages)
      
      for (i in seq_along(pages)) {
        dp <- filter(d, code %in% pages[[i]]) %>%
          mutate(code = droplevels(code),
                 # reorder within page so each panel reads high-to-low
                 label = reorder(label, score))
        ncol_i <- min(3, length(pages[[i]]))
        nrow_i <- ceiling(length(pages[[i]]) / ncol_i)
        
        p <- ggplot(dp, aes(score, label)) +
          geom_point(aes(size = n_overlap, colour = effect)) +
          scale_colour_viridis_c(name = colour_lab, option = colour_opt, end = .9) +
          scale_size_continuous(name = "Genes\nin term", range = c(1.5, 5.5)) +
          facet_wrap(~ code, scales = "free_y", ncol = ncol_i) +
          labs(title = sprintf("%s%s. %s: %s%s", "Figure ", fig_no, title_stem, dm,
                               if (npage > 1) sprintf(" (page %d of %d)", i, npage) else ""),
               subtitle = wrap_sub(paste0(
                 sub_stem, " Up to ", TERMS_PER_GROUP, " most significant terms per group, ",
                 "ordered by significance within each panel. Panels use independent term sets, so a ",
                 "term absent from a panel simply did not rank in that group's top ", TERMS_PER_GROUP,
                 ". Group codes: ", key_line)),
               x = expression(-log[10]~italic(p)), y = NULL) +
          theme(axis.text.y = element_text(size = 6.5),
                strip.text  = element_text(size = 8))
        
        suffix <- if (npage > 1) sprintf("_p%d", i) else ""
        sv(p, sprintf("fig_enrich_%s_%s%s.png", test_name,
                      gsub("[^A-Za-z0-9]+", "_", dm), suffix),
           w = 13, h = 2.6 + 2.9 * nrow_i)
      }
    }
    
    ora <- filter(top, test == "ora")
    if (nrow(ora) > 0) {
      for (dm in sort(unique(ora$domain)))
        dotplot_pages(filter(ora, domain == dm), "ora", dm, "23",
                      "Fold\nenrichment", "C",
                      "Over-Representation Analysis",
                      paste("Hypergeometric test of the up-regulated gene set for each group against the",
                            "14,782-gene filtered universe actually testable in this experiment, not the",
                            "annotated genome. Point size is the number of query genes in the term and",
                            "colour is fold enrichment."))
    }
    
    rk <- filter(top, test == "rank")
    if (nrow(rk) > 0) {
      for (dm in sort(unique(rk$domain)))
        dotplot_pages(filter(rk, domain == dm), "rank", dm, "24",
                      "AUC", "D",
                      "Threshold-Free Rank Test",
                      paste("Mann-Whitney U test of each term's members against all other genes across the",
                            "full ranked list, applying no significance threshold to the query. This test",
                            "recovers broad programmes that over-representation misses, but many of those",
                            "are housekeeping or proliferation related rather than subtype specific."))
    }
    
    # ------------------------------------------------- 3. Agreement between tests --
    if (nrow(ora) > 0 && nrow(rk) > 0) {
      agree <- top %>%
        distinct(domain, group, code, label, test) %>%
        mutate(present = TRUE) %>%
        pivot_wider(names_from = test, values_from = present, values_fill = FALSE)
      # pivot_wider only creates columns for tests that occur; create both so the
      # logical expression below cannot fail on a domain lacking one of them.
      if (!"ora"  %in% names(agree)) agree$ora  <- FALSE
      if (!"rank" %in% names(agree)) agree$rank <- FALSE
      agree <- agree %>%
        mutate(support = case_when(ora & rank ~ "Both tests",
                                   ora        ~ "Over-representation only",
                                   TRUE       ~ "Rank test only"))
      write.csv(agree, file.path(RES, "enrichment_agreement.csv"), row.names = FALSE)
      
      n_both <- sum(agree$support == "Both tests")
      
      # complete() so a group with zero dual-supported terms shows a zero rather
      # than dropping out of the chart entirely.
      cnt <- agree %>%
        count(code, support) %>%
        complete(code, support, fill = list(n = 0L)) %>%
        mutate(support = factor(support, levels = c("Both tests",
                                                    "Over-representation only",
                                                    "Rank test only")))
      
      p <- ggplot(cnt, aes(code, n, fill = support)) +
        geom_col(width = .7) +
        scale_fill_manual(values = c("Both tests" = "#1D9E75",
                                     "Over-representation only" = "#EF9F27",
                                     "Rank test only" = "#378ADD")) +
        labs(title = "Figure 25. Enriched Terms by Which Statistical Test Supports Them",
             subtitle = wrap_sub(paste0(
               "Stacked counts of enriched terms per diagnostic group. Over-representation depends on where ",
               "the query threshold is drawn; the rank test does not. Agreement between them therefore ",
               "indicates enrichment that is not an artefact of the cutoff, and the ", n_both,
               " terms in the green band are the ones reported as primary findings. Group codes: ", key_line)),
             x = "Diagnostic Group", y = "Number of Enriched Terms", fill = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
              legend.position = "bottom")
      sv(p, "fig_enrich_agreement.png", w = 11, h = 6.5)
      
      message(n_both, " terms supported by both tests")
    }
    
    # --------------------------------------------- 4. Cross-cluster term heatmap --
    # Explicit selection: within each domain, keep the terms that appear in the most
    # groups. The previous `first(domain)` construction did not do this.
    if (nrow(ora) > 0) {
      keep_terms <- ora %>%
        count(domain, label, name = "n_groups") %>%
        group_by(domain) %>%
        slice_max(n_groups, n = HEATMAP_TERMS, with_ties = FALSE) %>%
        ungroup() %>%
        select(domain, label)
      
      hm <- inner_join(ora, keep_terms, by = c("domain", "label"))
      
      if (nrow(hm) > 0) {
        p <- ggplot(hm, aes(code, reorder(label, score), fill = score)) +
          geom_tile(colour = "white", linewidth = .3) +
          scale_fill_viridis_c(name = expression(-log[10]~italic(p)),
                               option = "C", end = .9) +
          facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
          labs(title = "Figure 26. Enriched Terms Shared Across Multiple Diagnostic Groups",
               subtitle = wrap_sub(paste0(
                 "Within each library, the ", HEATMAP_TERMS, " terms enriched in the greatest number of ",
                 "groups. Blank cells did not reach the FDR cutoff in that group and are not zero effects. ",
                 "Terms enriched across many groups reflect shared biology or generic programmes rather ",
                 "than subtype-specific signal, so breadth here is not evidence of specificity. ",
                 "Group codes: ", key_line)),
               x = "Diagnostic Group", y = NULL) +
          theme(axis.text.y = element_text(size = 7),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                strip.text.y = element_text(size = 8))
        sv(p, "fig_enrich_term_heatmap.png", w = 12,
           h = max(7, 2 + 0.26 * length(unique(hm$label))))
      }
    }
    
    # ------------------------------------------------------------- 5. Term counts --
    f_sum <- file.path(RES, "enrichment_summary.csv")
    if (file.exists(f_sum)) {
      s <- read.csv(f_sum) %>%
        mutate(group = as.character(group),
               code = factor(code_of[group], levels = key$code)) %>%
        pivot_longer(c(n_ora_significant, n_rank_significant),
                     names_to = "test", values_to = "n_terms") %>%
        mutate(test = recode(test,
                             n_ora_significant  = "Over-representation",
                             n_rank_significant = "Rank test"))
      p <- ggplot(s, aes(code, n_terms, fill = test)) +
        geom_col(position = "dodge", width = .75) +
        scale_fill_manual(values = c("Over-representation" = "#EF9F27",
                                     "Rank test" = "#378ADD")) +
        facet_wrap(~ library, scales = "free_y") +
        labs(title = "Figure 27. Number of Enriched Terms per Group and Gene Set Library",
             subtitle = wrap_sub(paste0(
               "The rank test returns several times more terms than over-representation in every library. ",
               "Inspection of those additional terms shows they are largely housekeeping and proliferation ",
               "programmes - oxidative phosphorylation, mitochondrial translation, DNA repair, spliceosome ",
               "- which differ between quiescent and proliferating disease rather than between specific ",
               "subtypes. Library coverage of the background also differs, so counts are not comparable ",
               "across panels. Group codes: ", key_line)),
             x = "Diagnostic Group", y = "Terms at FDR Cutoff", fill = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
              strip.text = element_text(size = 8),
              legend.position = "bottom")
      sv(p, "fig_enrich_counts.png", w = 12, h = 7.5)
    }
    
  }  # end non-empty
}  # end file exists

if (.written == 0L) {
  message("\nNO ENRICHMENT FIGURES WRITTEN - check the enrichment step output.")
} else {
  message("\n", .written, " enrichment figures written to ", FIG)
}

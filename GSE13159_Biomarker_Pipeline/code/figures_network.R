#!/usr/bin/env Rscript
# =============================================================================
# Gene correlation network figures (Figures 20-22) - pipeline step 14
#
#   install.packages(c("igraph","ggplot2","dplyr","tidyr","RColorBrewer"))
#   Rscript figures_network.R
#
# FOUR ISSUES THIS FILE HAS TO HANDLE:
#
# 1. NEVER CALL quit() IN A SOURCED SCRIPT. The previous version called
#    quit(save="no") when network output was absent, which terminates the whole
#    R session - losing everything in the workspace when this file is source()d
#    from RStudio rather than run through Rscript. Missing input now skips the
#    section and returns normally.
#
# 2. PALETTE INDEXING BY MODULE ID IS UNSAFE. pal[module] assumes module IDs are
#    exactly 1..n. Any vertex absent from the nodes table, or any gap in the
#    numbering, silently yields NA and igraph then draws an uncoloured vertex.
#    Modules are mapped through a factor so every vertex gets a defined colour.
#
# 3. EDGE ATTRIBUTES MUST TRAVEL WITH THE GRAPH. Passing edge.color from a
#    separate data frame relies on igraph preserving row order exactly. It
#    currently does, but the assumption is invisible and breaks silently if the
#    edge list is ever reordered or de-duplicated. Sign is attached as a proper
#    edge attribute and read back off the graph.
#
# 4. DIVISION BY max(deg). An isolated subgraph gives max(deg) == 0 and vertex
#    sizes become NaN, which igraph renders as invisible points.
#
# INTERPRETIVE NOTE carried in the captions: on this cohort the global network
# is dominated by between-subtype contrast rather than gene-level regulation,
# so the within-cluster network is the interpretable one. Two hub genes of the
# small residual modules are sex-linked and reflect donor sex, not disease.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(RColorBrewer)
})
has_igraph <- requireNamespace("igraph", quietly = TRUE)

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
                  plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9, colour = "grey25")))

.written <- 0L
.missing <- character(0)

rd <- function(f) {
  p <- file.path(RES, f)
  if (!file.exists(p)) { .missing <<- c(.missing, f); return(NULL) }
  read.csv(p)
}
noted <- function(f) { .written <<- .written + 1L; message("wrote ", f) }
sv <- function(p, f, w = 9, h = 5) {
  ggsave(file.path(FIG, f), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  noted(f)
}
wrap_sub <- function(x, width = 110) paste(strwrap(x, width = width), collapse = "\n")

nodes <- rd("08_network_nodes.csv")

# Everything below is inside this guard rather than after a quit() call, so a
# missing input skips the section instead of destroying the caller's session.
if (is.null(nodes)) {
  message("no network output found - run s08_gene_network.py to produce ",
          "08_network_nodes.csv. Skipping the network figures.")
} else {
  
  n_gene <- nrow(nodes)
  
  # ---------------------------------------------------- 1. degree distribution
  if (all(c("degree_global", "degree_within") %in% names(nodes))) {
    tot_g <- sum(nodes$degree_global) / 2
    tot_w <- sum(nodes$degree_within) / 2
    pct   <- if (tot_g > 0) 100 * tot_w / tot_g else NA_real_
    
    d <- nodes %>%
      pivot_longer(c(degree_global, degree_within),
                   names_to = "network", values_to = "degree") %>%
      mutate(network = recode(network,
                              degree_global = "Global (all samples pooled)",
                              degree_within = "Within-cluster"))
    p <- ggplot(d, aes(degree, fill = network)) +
      geom_histogram(bins = 40, position = "identity", alpha = .6) +
      scale_fill_manual(values = c("Global (all samples pooled)" = "#C44E52",
                                   "Within-cluster" = "#4C72B0")) +
      labs(title = "Figure 20. Node Degree Distribution Before and After Removing Between-Subtype Contrast",
           subtitle = wrap_sub(sprintf(paste(
             "Degree is the number of retained correlation partners per gene. Computing correlations",
             "across the pooled cohort (red) yields %s edges; recomputing them within clusters and",
             "aggregating (blue) yields %s, or %.1f%% of that total. The difference is not noise removal",
             "but the removal of between-subtype contrast: when transcriptionally distant subtypes are",
             "pooled, thousands of genes shift together and pairwise correlation measures group identity",
             "rather than gene-level regulation."),
             format(tot_g, big.mark = ","), format(tot_w, big.mark = ","), pct)),
           x = "Node Degree", y = "Number of Genes", fill = NULL) +
      theme(legend.position = "bottom")
    sv(p, "fig_net_degree.png", w = 9.5, h = 5.5)
  }
  
  # ---------------------------------------------------------- 2. module sizes
  if ("module" %in% names(nodes)) {
    ms <- count(nodes, module)
    biggest <- ms$n[which.max(ms$n)]
    p <- ggplot(ms, aes(factor(module), n)) +
      geom_col(fill = "#4C72B0", width = .7) +
      geom_text(aes(label = n), vjust = -0.4, size = 3) +
      labs(title = "Figure 21. Size of Each Detected Correlation Module in the Global Network",
           subtitle = wrap_sub(sprintf(paste(
             "Module membership from community detection on the global correlation network. One module",
             "contains %d of %d genes (%.0f%%), which is not a co-expression module in any useful sense:",
             "it is the signature of a near-complete graph produced by pooling distinct subtypes. The",
             "remaining modules are small and are interpreted with caution, since two of their hub genes",
             "are sex-linked and therefore reflect donor sex rather than disease biology."),
             biggest, n_gene, 100 * biggest / n_gene)),
           x = "Module", y = "Number of Genes")
    sv(p, "fig_net_modules.png", w = 8.5, h = 5.5)
  }
  
  # ------------------------------------------------------ 3. network layouts
  panel <- c(global = "22A", within = "22B")
  net_desc <- c(
    global = paste("Correlations computed across all samples. Interpret with caution: the giant",
                   "component reflects subtype identity rather than gene-level regulation."),
    within = paste("Correlations computed within clusters and aggregated, which removes",
                   "between-subtype contrast. This is the interpretable network."))
  
  for (which_net in c("global", "within")) {
    e <- rd(paste0("08_network_edges_", which_net, ".csv"))
    if (is.null(e) || nrow(e) == 0) next
    if (!has_igraph) {
      message("igraph is not installed - skipping the layout plots. ",
              "install.packages(\"igraph\") to enable them.")
      break
    }
    
    n_edge_all <- nrow(e)
    shown_note <- ""
    if (n_edge_all > 4000) {
      e <- e[order(-abs(e$r)), ][1:4000, ]
      shown_note <- sprintf(paste(" For legibility only the 4,000 strongest of %s edges are drawn,",
                                  "so absolute connectivity is understated here; see Figure 20 for",
                                  "the full degree distribution."), format(n_edge_all, big.mark = ","))
      message(which_net, ": showing the 4000 strongest of ", n_edge_all, " edges")
    }
    
    g <- igraph::graph_from_data_frame(e[, c("gene_a", "gene_b")], directed = FALSE)
    
    # Attach sign as a true edge attribute so colour cannot drift out of step
    # with the edge list if it is ever reordered.
    if ("sign" %in% names(e)) igraph::E(g)$sign <- as.character(e$sign)
    else                      igraph::E(g)$sign <- "positive"
    ecol <- ifelse(igraph::E(g)$sign == "positive", "#C4525288", "#4C72B088")
    
    # Map modules through a factor: indexing a palette by raw module ID gives
    # NA for any vertex missing from the nodes table or any gap in numbering.
    mod_lookup <- setNames(nodes$module, nodes$gene)
    vmod <- mod_lookup[igraph::V(g)$name]
    vmod[is.na(vmod)] <- -1L                     # explicit "unassigned" level
    vfac <- factor(vmod)
    pal  <- colorRampPalette(brewer.pal(8, "Set2"))(nlevels(vfac))
    vcol <- pal[as.integer(vfac)]
    
    deg  <- igraph::degree(g)
    dmax <- max(deg, 1)                          # guard against max(deg) == 0
    vsize <- 2 + 6 * deg / dmax
    lab_cut <- if (length(deg)) quantile(deg, .97) else Inf
    
    png(file.path(FIG, paste0("fig_net_layout_", which_net, ".png")),
        width = 2000, height = 2150, res = 200)
    op <- par(mar = c(1, 1, 7, 1))
    set.seed(42)
    plot(g, layout = igraph::layout_with_fr(g),
         vertex.color = vcol, vertex.frame.color = NA,
         vertex.size = vsize,
         vertex.label = ifelse(deg >= lab_cut, igraph::V(g)$name, NA),
         vertex.label.cex = .6, vertex.label.color = "black",
         edge.color = ecol, edge.width = .4,
         main = "")
    title(main = paste0("Figure ", panel[[which_net]],
                        ". Gene Correlation Network (", which_net,
                        "), Vertices Coloured by Module"),
          cex.main = 1.15, font.main = 2, line = 5)
    mtext(paste(strwrap(paste0(
      net_desc[[which_net]],
      " Vertex area scales with degree and the highest-degree 3% of genes are labelled.",
      " Edge colour indicates the sign of the correlation (red positive, blue negative).",
      shown_note), width = 118), collapse = "\n"),
      side = 3, line = 0.4, cex = 0.72, col = "grey25", adj = 0)
    par(op)
    dev.off()
    noted(paste0("fig_net_layout_", which_net, ".png"))
  }
}

# ------------------------------------------------------------------ summary --
if (length(.missing)) {
  message("\nskipped, inputs not found in ", RES, ":\n  ",
          paste(unique(.missing), collapse = "\n  "))
}
if (.written == 0L) {
  message("\nNO NETWORK FIGURES WRITTEN - check that s08_gene_network.py has run.")
} else {
  message("\n", .written, " network figures written to ", FIG)
}

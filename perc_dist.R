perc_dist <- function(condor, anno){
perc_tab <- table(condor$anno$cell_anno[,anno])/length(condor$anno$cell_anno[,anno])*100
return(perc_tab)}


# =============================================================================
# condor_cluster_composition
#
# Plots cell composition per cluster, grouped by an annotation variable.
# Error bars = SD across replicates (samples) within each group.
# Outlier points = individual sample values overlaid.
#
# Arguments:
#   condor        : condor object
#   anno_col      : annotation column to group by
#                   e.g. "condition", "timepoint"
#                   must be present in condor$anno$cell or joinable via sample_id
#   metric        : "percentage" or "count"
#   sid_col       : sample ID column name
#   min_cells     : clusters with fewer total cells than this are dropped
#   errorbar_fun  : "sd" or "se" (standard error)
#   palette       : named or unnamed character vector of colours; NULL = default
# =============================================================================

condor_cluster_composition <- function(
    condor,
    cluster_slot,
    clust_var = NULL,
    cluster_col = NULL,
    anno_col,
    metric       = c("percentage", "count"),
    sid_col      = "sample_id",
    min_cells    = 50,
    errorbar_fun = c("sd", "se"),
    palette      = NULL
) {
  
  metric       <- match.arg(metric)
  errorbar_fun <- match.arg(errorbar_fun)
  
  # ── dependencies ────────────────────────────────────────────────────────────
  required <- c("ggplot2","dplyr","tidyr","ggrepel")
  missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0)
    stop("Install missing packages: ", paste(missing, collapse = ", "))
  
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(tidyr)
  })  
  # ── extract and join data ────────────────────────────────────────────────────
  if (is.null(condor$clustering))
    stop("condor$clustering is NULL — run clustering first.")
  
  if (!is.list(condor$clustering))
    stop("condor$clustering is not a list. Structure may differ from expected.")
  
  # show available resolutions if cluster_col not found
  available <- names(condor$clustering)
  
  if (!cluster_slot %in% available)
    stop("'", cluster_slot, "' not found in condor$clustering.\n",
         "Available resolutions: ", paste(available, collapse = ", "))
  
  clust_df <- condor$clustering[[cluster_slot]]
  colnames(clust_df)[1:2] <- c("cluster","Description")
  clust_df[,"cluster"] <- condor$clustering[[cluster_slot]][,clust_var]
  
  # ── join annotation ──────────────────────────────────────────────────────────
  if (!anno_col %in% colnames(clust_df)) {
    
    # check condor$anno$cell_anno first
    cell_anno <- condor$anno$cell_anno
    
    if (is.null(cell_anno))
      stop("condor$anno$cell_anno is NULL.")
    
    if (!anno_col %in% colnames(cell_anno))
      stop("'", anno_col, "' not found in condor$anno$cell_anno.\n",
           "Available columns: ", paste(colnames(cell_anno), collapse = ", "))
    
    if (!sid_col %in% colnames(cell_anno))
      stop("'", sid_col, "' not found in condor$anno$cell_anno.")
    
    anno_lookup <- cell_anno %>%
      select(all_of(c(sample=sid_col, group=anno_col)))
    
    clust_df <- merge(clust_df, anno_lookup, by = "row.names")
  }
  
  if (!is.data.frame(clust_df))
    stop("condor$clustering[['", cluster_slot, "']] is not a dataframe.")
  # ── warn if error bars will be uninformative ─────────────────────────────────
  n_replicates <- clust_df %>%
    distinct(group, sample) %>%
    count(group) %>%
    pull(n)
  
  if (any(n_replicates < 2))
    warning("Some groups have only 1 sample — error bars will be NA for those groups. ",
            "Consider using a grouping variable with multiple samples per level.")
  
  # ── drop small clusters ──────────────────────────────────────────────────────
  cluster_totals <- clust_df %>%
    count(cluster, name = "total_cells")
  
  small_clusters <- cluster_totals %>%
    filter(total_cells < min_cells) %>%
    pull(cluster)
  
  if (length(small_clusters) > 0)
    message("Dropping ", length(small_clusters), " cluster(s) with < ",
            min_cells, " cells: ", paste(small_clusters, collapse = ", "))
  
  clust_df <- clust_df %>%
    filter(!cluster %in% small_clusters)
  
  # ── compute per-sample counts and percentages ────────────────────────────────
  # count cells per sample × cluster
  per_sample <- clust_df %>%
    count(sample, group, cluster, name = "n_cells")
  
  # complete: ensure every sample × cluster combination exists (zero if absent)
  per_sample <- per_sample %>%
    complete(nesting(sample, group), cluster,
             fill = list(n_cells = 0))
  
  # compute percentage within each sample (of total cells in that sample)
  per_sample <- per_sample %>%
    group_by(sample) %>%
    mutate(
      total_in_sample = sum(n_cells),
      percentage      = ifelse(total_in_sample > 0,
                               n_cells / total_in_sample * 100,
                               0)
    ) %>%
    ungroup()
  
  # choose the value column based on metric
  per_sample <- per_sample %>%
    mutate(value = if (metric == "percentage") percentage else n_cells)
  
  # ── summary stats per group × cluster ───────────────────────────────────────
  summary_df <- per_sample %>%
    group_by(group, cluster) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val   = sd(value,   na.rm = TRUE),
      n_rep    = n(),
      se_val   = sd_val / sqrt(n_rep),
      .groups  = "drop"
    ) %>%
    mutate(
      err = if (errorbar_fun == "sd") sd_val else se_val,
      ymin = pmax(mean_val - err, 0),   # floor at 0
      ymax = mean_val + err
    )
  
  # ── axis label ───────────────────────────────────────────────────────────────
  y_label <- if (metric == "percentage") {
    paste0("% of cells per sample (± ", errorbar_fun, ")")
  } else {
    paste0("Cell count per sample (± ", errorbar_fun, ")")
  }
  
  # ── colour palette ───────────────────────────────────────────────────────────
  groups <- sort(unique(summary_df$group))
  n_groups <- length(groups)
  
  if (is.null(palette)) {
    default_pal <- c("#185FA5","#D85A30","#1D9E75","#7F77DD",
                     "#EF9F27","#A32D2D","#5DCAA5","#888780",
                     "#3CBFE0","#F4A261","#2D6A4F","#C77DFF")
    palette <- setNames(default_pal[seq_len(n_groups)], groups)
  } else if (is.null(names(palette))) {
    palette <- setNames(rep_len(palette, n_groups), groups)
  }
  
  # ── plot ─────────────────────────────────────────────────────────────────────
  p <- ggplot(summary_df,
              aes(x     = cluster,
                  y     = mean_val,
                  fill  = group,
                  group = group)) +
    
    # bars
    geom_bar(stat     = "identity",
             position = position_dodge(width = 0.8),
             width    = 0.7,
             alpha    = 0.85) +
    
    # error bars
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      position = position_dodge(width = 0.8),
      width    = 0.25,
      linewidth = 0.6,
      colour   = "grey30"
    ) +
    
    # individual sample points (outliers visible)
    geom_point(
      data     = per_sample,
      aes(x     = cluster,
          y     = value,
          group = group,
          colour = group),
      position = position_dodge(width = 0.8),
      size     = 1.8,
      shape    = 21,
      fill     = "white",
      stroke   = 0.8,
      show.legend = FALSE
    ) +
    
    scale_fill_manual(values = palette, name = anno_col) +
    scale_colour_manual(values = palette, name = anno_col) +
    
    labs(
      title    = paste("Cluster composition by", anno_col),
      subtitle = paste0(
        metric, "  |  error bars = ", errorbar_fun,
        "  |  points = individual samples",
        "  |  clustering = ", cluster_slot
      ),
      x = "Cluster",
      y = y_label
    ) +
    
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 9, colour = "grey50"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "top"
    )
  
  print(p)
}

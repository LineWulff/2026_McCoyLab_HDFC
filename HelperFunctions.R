#' Helper script for R script CyCondor pipeline for high dimensional flow cytometry data
#' Includes all helper functions called from main script
#' Author: Line Wulff
#' Date (created): 26-06-09

# =============================================================================
# perc_dist
#
# Returns the composition of cells from your specificed annotation
#
# Arguments:
#   condor        : condor object
#   anno      : annotation column to group by
#                   e.g. "condition", "timepoint"
#                   must be present in condor$anno$cell or joinable via sample_id
# =============================================================================
perc_dist <- function(condor, anno){
  perc_tab <- table(condor$anno$cell_anno[,anno])/length(condor$anno$cell_anno[,anno])*100
  return(perc_tab)}

# =============================================================================
# find_high_cor_pairs
#
# Returns marker pairs with correlation above specified threshold for inspection
#
# Arguments:
#   cor_mat       : correlation matrix of markers
#   threshold      : threshold for correlation pairs to be returned, set automatically
#                   to 0.5, but down to 0.2 is a good idea for proper QC check
# =============================================================================
find_high_cor_pairs <- function(cor_mat, threshold = 0.5) {
  
  # Basic checks
  if (!is.matrix(cor_mat)) stop("Input must be a matrix")
  if (is.null(rownames(cor_mat)) || is.null(colnames(cor_mat))) {
    stop("Matrix must have row and column names")
  }
  
  # Get upper triangle indices (avoid duplicates + diagonal)
  upper_idx <- upper.tri(cor_mat, diag = FALSE)
  
  # Extract values
  vals <- cor_mat[upper_idx]
  
  # Get row/col names for those positions
  row_names <- rownames(cor_mat)[row(cor_mat)[upper_idx]]
  col_names <- colnames(cor_mat)[col(cor_mat)[upper_idx]]
  
  # Filter by threshold (exclude 1s automatically since diag removed)
  keep <- vals > threshold
  
  # Return as data frame
  data.frame(
    feature1 = row_names[keep],
    feature2 = col_names[keep],
    correlation = vals[keep],
    row.names = NULL
  )
}

# =============================================================================
# plot_pseudobulkPCA
#
# Makes a PCA plot coloured/shaped points by up to two sample annotations
#
#   condor        : condor object
#   extra_anno      : annotation column to group by
#                   e.g. "condition", "timepoint"
#                   must be present in condor$anno$cell or joinable via sample_id
#                   extra_anno accepts up to two arguments, if you don't want to color 
#                   or shape by anything set extra_anno = NULL
#   marker_cols     : which markers to use for the PCA
#   AnnoTable      : Your sample annotation table
# =============================================================================
plot_pseudobulkPCA <- function(condor, 
                               extra_anno = NULL, 
                               marker_cols = colnames(condor$expr$orig), 
                               AnnoTable=AnnoTable){
  med_mat <- cbind(condor$expr$orig, sample_id=condor$anno$cell_anno$sample_id) %>%
    group_by(sample_id) %>%
    summarise(across(all_of(marker_cols), median))
  
  pca <- prcomp(med_mat[, marker_cols], scale.=TRUE)
  pca_df <- as.data.frame(pca$x[, 1:2])
  pca_df$sample    <- med_mat$sample_id
  if (length(extra_anno)>0){
    for (i in seq(1,length(extra_anno))){
      ann <- extra_anno[i]
      pca_df[,ann] <- AnnoTable[pca_df$sample, ann]}} else {ann1 = NULL: ann2 = NULL}
  print(pca_df)
  point_aes <- switch(as.character(length(extra_anno)),
                      "0" = aes(x = PC1, y = PC2, label = sample),
                      "1" = aes(x = PC1, y = PC2,label = sample,
                                colour = .data[[extra_anno[1]]]),
                      "2" = aes(x = PC1, y = PC2, label = sample,
                                colour = .data[[extra_anno[1]]],
                                shape  = .data[[extra_anno[2]]]))
  
  pct <- round(summary(pca)$importance[2, 1:2]*100, 1)
  pca_pseu_plot <- ggplot(pca_df, point_aes) +
    geom_point(size=4) +
    ggrepel::geom_text_repel(size=3) +
    labs(title="PCA of sample medians — pre-UMAP batch check",
         x=paste0("PC1 (",pct[1],"%)"),
         y=paste0("PC2 (",pct[2],"%)")) +
    theme_minimal()
  return(pca_pseu_plot)
}

# =============================================================================
# geom_cluster_labels
#
# adds cluster labels to your ggplot based UMAP
#
# Arguments:
#   data        : the same dataframe you used to plot your UMAP
#   cluster_col     : column name of cluster identities ou want to colour by
#   x_col           : your x axis, UMAP1 if bot specified, has to be the same as in main plot
#   y_col           : your y axis, UMAP2 if not specified, has to be the same as in main plot
#   label_size      : size of labels
#   seed            : for reproducibility (helps decides exact labe position)
# =============================================================================
geom_cluster_labels <- function(
    data,                      # the same df used in ggplot()
    cluster_col,               # string: column name with cluster identities
    x_col       = "UMAP1",    # string: column name for x axis
    y_col       = "UMAP2",    # string: column name for y axis
    label_size  = 4,
    seed        = 42
) {
  
  # ── compute centroids ──────────────────────────────────────────────────────
  centroids <- data %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(
      UMAP1 = median(.data[[x_col]], na.rm = TRUE),
      UMAP2 = median(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    )
  colnames(centroids)[1] <- "label"
  
  # ── build layer list ───────────────────────────────────────────────────────
  layers <- list()
  # # simplest: fixed text at centroid, no box
  layers[[1]] <- geom_text(
    data        = centroids,
    mapping     = aes(x = get(x_col),
                      y = get(y_col)),
    label       = centroids$label,   # passed outside aes() entirely
    size        = label_size,
    fontface    = "bold",
    colour      = "black",
    inherit.aes = FALSE
  )
  
  layers
}

# =============================================================================
# runPhenograph_issue35
#
# Some systems have errors w. phenograph, if you get error message:
# Error in if (any(i < 0L)) { : missing value where TRUE/FALSE needed
# This function is written by the cycondor developpers and solves this issue
#
# Arguments:
#   see runPhenograph - same arguments needed
# =============================================================================
runPhenograph_issue35 <- function (fcd,
                                   input_type,
                                   data_slot,
                                   k,
                                   seed = 91,
                                   prefix = NULL,
                                   nPC = ncol(fcd$pca[[data_slot]]),
                                   markers = colnames(fcd$expr[[data_slot]]),
                                   discard = FALSE){
  set.seed(seed)
  
  # see if selected markers are present in condor_object , if input_type = "expr"
  if (input_type == "expr"){
    for (single in markers){
      if (!single %in% colnames(fcd[["expr"]][["orig"]])){
        stop(paste("ERROR:",single, "not found in expr markers."))
      }
    }
    
    # define markers to use
    if (discard == FALSE){              # (discard == F -> keep specified markers (default = all))
      
      phe_markers <- markers
    }
    
    else if (discard == TRUE) {       # (discard == T -> discard specified markers, error if no markers are specified)
      
      if (length(markers) == length(colnames(fcd$expr[[data_slot]]))){
        
        stop("ERROR: No markers specified. Specify markers to be removed or set 'discard = F'.")
      }
      
      else {
        
        phe_markers <- setdiff(colnames(fcd$expr[[data_slot]]), markers)
        
      }
    }
    
    #define fcd subset for calculation
    data1 <- fcd$expr[[data_slot]][, colnames(fcd$expr[[data_slot]]) %in% phe_markers, drop = F]
    
    
    
  }
  if (input_type == "pca"){
    
    #define fcd subset for calculations and get used markers of PCA analysis
    
    data1 <- fcd$pca[[data_slot]][,1:nPC]
    phe_markers <- used_markers(fcd,  input_type = "pca", data_slot = data_slot, mute = T)
  }
  
  
  
  
  Rphenograph_out <- Rphenograph::Rphenograph(data1, k = k)
  Rphenograph_out <- as.matrix(igraph::membership(Rphenograph_out[[2]]))
  Rphenograph_out <- as.data.frame(matrix(ncol = 1,
                                          data = Rphenograph_out,
                                          dimnames = list(rownames(fcd$expr$orig), "Phenograph")))
  Rphenograph_out$Phenograph <- as.factor(Rphenograph_out$Phenograph)
  Rphenograph_out$Description <- paste(input_type, "_", data_slot, "_k", k, sep = "")
  
  phe_name <- paste("phenograph", sub("^_", "", paste(prefix, input_type, data_slot, "k", k, sep = "_")), sep = "_")
  
  if (input_type == "pca"){
    if (nPC < ncol(fcd[[input_type]][[data_slot]])) {
      
      suffix <- paste0("top", nPC)
      
      phe_name <- paste(phe_name, suffix, sep = "_")
    }
  }
  
  fcd[["clustering"]][[phe_name]] <- Rphenograph_out
  
  #save used markers in "extras"-slot
  fcd[["extras"]][["markers"]][[paste(phe_name,"markers", sep = "_")]] <- phe_markers
  
  
  return(fcd)
}

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

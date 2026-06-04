geom_cluster_labels <- function(
    data,                      # the same df used in ggplot()
    cluster_col,               # string: column name with cluster identities
    x_col       = "UMAP1",    # string: column name for x axis
    y_col       = "UMAP2",    # string: column name for y axis
    label_size  = 4,
    # label_box   = TRUE,        # TRUE = geom_label style, FALSE = geom_text
    # repel       = TRUE,        # TRUE = ggrepel, FALSE = fixed at centroid
    # palette     = NULL,        # named vector of colours; NULL = inherit from plot
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
  # # ── compute centroids ──────────────────────────────────────────────────────
  # centroids <- as.data.frame(
  #   do.call(rbind, lapply(unique(data[[cluster_col]]), function(cl) {
  #     sub <- data[data[[cluster_col]] == cl, ]
  #     data.frame(
  #       label = as.character(cl),
  #       x     = median(sub[[x_col]], na.rm = TRUE),
  #       y     = median(sub[[y_col]], na.rm = TRUE),
  #       stringsAsFactors = FALSE
  #     )
  #   }))
  # )
  
  # ── resolve fill colours ───────────────────────────────────────────────────
  # if palette supplied, map cluster → colour for filled labels
  # if NULL, fall back to dark text with white shadow (works on any bg)
  # use_fill <- !is.null(palette)
  # 
  # if (use_fill) {
  #   if (is.null(names(palette))) {
  #     # unnamed vector — map by factor order
  #     lvls <- levels(factor(data[[cluster_col]]))
  #     palette <- setNames(rep_len(palette, length(lvls)), lvls)
  #   }
  #   centroids$fill_col <- palette[as.character(centroids$label)]
  # }
  
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

# =============================================================================
# cyCONDOR post-load QC wrapper
# Runs all 8 QC checks and saves an HTML report + individual PDF plots
#
# Usage:
#   source("condor_qc.R")
#   run_condor_qc(
#     condor    = your_condor_object,
#     metadata  = your_metadata_df,      # rownames = sample IDs
#     out_dir   = "qc_output",
#     panel     = your_panel_df,         # optional: data.frame with $marker col
#     cofactor  = 150                    # arcsinh cofactor used at read-in
#   )
#
# metadata must have at minimum:
#   - rownames matching sample_id values in condor$expr$orig
#   - a "condition" column
#   - optionally: "sex", "batch", "timepoint"
#
# Output:
#   qc_output/
#     condor_qc_report.html   <- self-contained HTML report
#     plots/                  <- individual PNG files for each check
# =============================================================================

run_condor_qc <- function(condor,
                          metadata,
                          out_dir   = "qc_output",
                          panel     = NULL,
                          cofactor  = 150,
                          n_cells_sample = 5000,
                          seed = 42) {
  
  # ── dependencies ────────────────────────────────────────────────────────────
  required_pkgs <- c("ggplot2","dplyr","tidyr","ggridges",
                     "corrplot","ggrepel","patchwork","rmarkdown","knitr")
  missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly=TRUE)]
  if (length(missing) > 0) {
    stop("Please install missing packages:\n  install.packages(c(",
         paste0('"', missing, '"', collapse=", "), "))")
  }
  
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(tidyr)
    library(ggridges); library(corrplot); library(ggrepel)
    library(patchwork)
  })
  
  # ── setup ───────────────────────────────────────────────────────────────────
  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  set.seed(seed)
  theme_qc <- theme_minimal(base_size = 11) +
    theme(plot.title   = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey50"),
          strip.text   = element_text(size = 9))
  
  msgs   <- list()   # collects pass/warn/fail messages per check
  plots  <- list()   # collects ggplot objects for report
  
  save_plot <- function(p, name, w = 8, h = 5) {
    path <- file.path(plot_dir, paste0(name, ".png"))
    ggsave(path, p, width = w, height = h, dpi = 150, bg = "white")
    plots[[name]] <<- path
    invisible(path)
  }
  
  flag <- function(check, status, msg) {
    icon <- switch(status, pass = "PASS", warn = "WARN", fail = "FAIL")
    msgs[[check]] <<- list(status = status, icon = icon, msg = msg)
    cat(sprintf("[%s] %s: %s\n", icon, check, msg))
  }
  
  # ── helper: get expression matrix ───────────────────────────────────────────
  get_expr <- function() {
    if (!is.null(condor$expr$orig)) return(condor$expr$orig)
    if (!is.null(condor$expr$normalized)) return(condor$expr$normalized)
    stop("Cannot find condor$expr$orig or condor$expr$normalized")
  }
  
  expr      <- get_expr()
  sid_col   <- if ("sample_id" %in% colnames(expr)) "sample_id" else
    if ("expfcs_filename" %in% colnames(expr)) "expfcs_filename" else
      stop("Cannot identify sample ID column in condor expression slot")
  
  meta_cols  <- c(sid_col, "condition",
                  intersect(c("sex","batch","timepoint"), colnames(expr)))
  marker_cols <- setdiff(colnames(expr),
                         c(sid_col, "cell_id", "condition", "sex",
                           "batch", "timepoint", "sample_id",
                           "expfcs_filename"))
  # drop any residual scatter / time channels
  marker_cols <- marker_cols[!grepl("^(FSC|SSC|Time|time|Event)",
                                    marker_cols, ignore.case = TRUE)]
  
  cat(sprintf("\ncyCONDOR QC — %d cells, %d markers, %d samples\n\n",
              nrow(expr), length(marker_cols),
              length(unique(expr[[sid_col]]))))
  
  # ── 1. MARKER NAME & CHANNEL INTEGRITY ──────────────────────────────────────
  cat("── Check 1/8: Marker names & channel integrity\n")
  
  n_expected <- if (!is.null(panel)) nrow(panel) else NA
  n_found    <- length(marker_cols)
  na_counts  <- colSums(is.na(expr[, marker_cols]))
  n_na_mkrs  <- sum(na_counts > 0)
  scatter_in <- sum(grepl("^(FSC|SSC|Time)", colnames(condor$expr$orig %||% condor$expr$normalized), ignore.case=TRUE))
  
  status_1 <- if (n_na_mkrs > 0) "fail" else if (!is.na(n_expected) && abs(n_found - n_expected) > 2) "warn" else "pass"
  flag("marker_integrity", status_1,
       sprintf("%d markers found; %d with NAs; %s scatter channels in matrix",
               n_found,
               n_na_mkrs,
               if (scatter_in > 0) paste(scatter_in, "REMAINING") else "no"))
  
  # NA heatmap per sample × marker
  na_df <- expr %>%
    group_by(.data[[sid_col]]) %>%
    summarise(across(all_of(marker_cols), ~sum(is.na(.)), .names = "{.col}"),
              .groups = "drop") %>%
    pivot_longer(-1, names_to = "marker", values_to = "n_na")
  
  p1 <- ggplot(na_df, aes(x = marker, y = .data[[sid_col]], fill = n_na)) +
    geom_tile(colour = "white") +
    scale_fill_gradient(low = "white", high = "#A32D2D", name = "NA count") +
    labs(title  = "Check 1: Marker NA counts per sample",
         subtitle = "All white = no missing data",
         x = NULL, y = NULL) +
    theme_qc +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
  save_plot(p1, "01_marker_na_heatmap", w = 10, h = 4)
  
  # ── 2. METADATA JOIN INTEGRITY ───────────────────────────────────────────────
  cat("── Check 2/8: Metadata join integrity\n")
  
  expr <- expr %>%
    left_join(metadata %>% tibble::rownames_to_column(sid_col),
              by = sid_col, suffix = c("", ".meta"))
  
  if (!"condition" %in% colnames(expr)) {
    expr$condition <- NA_character_
  }
  
  n_na_cond   <- sum(is.na(expr$condition))
  unmatched   <- setdiff(unique(expr[[sid_col]]), rownames(metadata))
  status_2    <- if (n_na_cond > 0 || length(unmatched) > 0) "fail" else "pass"
  flag("metadata_join", status_2,
       sprintf("%d cells with NA condition; %d sample IDs not in metadata",
               n_na_cond, length(unmatched)))
  
  cross_df <- expr %>%
    count(.data[[sid_col]], condition) %>%
    mutate(condition = replace_na(as.character(condition), "NA — JOIN FAILED"))
  
  p2 <- ggplot(cross_df, aes(x = reorder(.data[[sid_col]], n),
                             y = n, fill = condition)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title    = "Check 2: Cells per sample after metadata join",
         subtitle = "Red 'NA' bars = metadata join failure",
         x = NULL, y = "Cell count") +
    scale_fill_manual(values = c("#185FA5","#1D9E75","#D85A30",
                                 "#A32D2D","#888780","#7F77DD"),
                      na.value = "#A32D2D") +
    theme_qc
  save_plot(p2, "02_metadata_join", w = 7, h = 4)
  
  # ── 3. TRANSFORMATION CORRECTNESS ───────────────────────────────────────────
  cat("── Check 3/8: Transformation correctness\n")
  
  sub_expr <- expr %>%
    slice_sample(n = min(n_cells_sample * length(unique(expr[[sid_col]])),
                         nrow(expr))) %>%
    select(all_of(c(sid_col, "condition", marker_cols)))
  
  long_expr <- sub_expr %>%
    pivot_longer(all_of(marker_cols), names_to = "marker", values_to = "intensity")
  
  range_stats <- long_expr %>%
    group_by(marker) %>%
    summarise(
      p05  = quantile(intensity, 0.05, na.rm = TRUE),
      p50  = median(intensity, na.rm = TRUE),
      p95  = quantile(intensity, 0.95, na.rm = TRUE),
      pmin = min(intensity, na.rm = TRUE),
      pmax = max(intensity, na.rm = TRUE),
      pct_below_neg2 = mean(intensity < -2, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  n_bad_range <- sum(range_stats$p95 < 0.5 | range_stats$pct_below_neg2 > 10)
  status_3 <- if (n_bad_range > 2) "fail" else if (n_bad_range > 0) "warn" else "pass"
  flag("transformation", status_3,
       sprintf("%d markers with suspicious range (p95 < 0.5 or >10%% below -2); cofactor used = %g",
               n_bad_range, cofactor))
  
  # Show up to 20 markers in ridgeline
  top_markers <- unique(long_expr$marker)[seq_len(min(20, length(marker_cols)))]
  p3 <- long_expr %>%
    filter(marker %in% top_markers) %>%
    ggplot(aes(x = intensity, y = marker, fill = marker)) +
    geom_density_ridges(alpha = 0.7, scale = 0.85, show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(title    = "Check 3: Transformed marker distributions",
         subtitle = sprintf("arcsinh cofactor = %g | each line = pooled cells", cofactor),
         x = "arcsinh intensity", y = NULL) +
    theme_qc
  save_plot(p3, "03_transformation_ridges", w = 8, h = max(5, length(top_markers) * 0.4 + 2))
  
  # ── 4. DOWNSAMPLING BALANCE ──────────────────────────────────────────────────
  cat("── Check 4/8: Downsampling balance\n")
  
  counts_df <- expr %>%
    count(.data[[sid_col]], condition) %>%
    mutate(pct = n / sum(n) * 100)
  
  min_pct    <- min(counts_df$pct)
  max_pct    <- max(counts_df$pct)
  ratio      <- max_pct / min_pct
  status_4   <- if (ratio > 4) "fail" else if (ratio > 2) "warn" else "pass"
  flag("downsampling", status_4,
       sprintf("min sample contribution %.1f%%, max %.1f%% (ratio %.1fx)",
               min_pct, max_pct, ratio))
  
  p4 <- ggplot(counts_df, aes(x = reorder(.data[[sid_col]], n),
                              y = n, fill = condition)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = mean(counts_df$n),
               linetype = "dashed", colour = "grey40") +
    coord_flip() +
    labs(title    = "Check 4: Cell counts per sample in condor object",
         subtitle = "Dashed line = mean; bars should be roughly equal",
         x = NULL, y = "Cell count") +
    theme_qc
  save_plot(p4, "04_downsampling_balance", w = 7, h = 4)
  
  # ── 5. PER-MARKER OUTLIERS ───────────────────────────────────────────────────
  cat("── Check 5/8: Per-marker outlier check\n")
  
  outlier_stats <- sapply(marker_cols, function(m) {
    x   <- expr[[m]]
    p95 <- quantile(x, 0.95, na.rm = TRUE)
    p99 <- quantile(x, 0.99, na.rm = TRUE)
    mx  <- max(x, na.rm = TRUE)
    ratio <- if (p95 > 0.1) p99 / p95 else NA_real_
    c(p95 = p95, p99 = p99, max = mx, ratio = ratio)
  })
  outlier_df <- as.data.frame(t(outlier_stats))
  outlier_df$marker <- rownames(outlier_df)
  outlier_df <- outlier_df[order(-outlier_df$ratio, na.last = TRUE), ]
  
  n_outlier_mkrs <- sum(outlier_df$ratio > 2, na.rm = TRUE)
  status_5 <- if (n_outlier_mkrs > 3) "warn" else "pass"
  flag("outliers", status_5,
       sprintf("%d markers with 99th/95th percentile ratio > 2 (winsorising recommended)",
               n_outlier_mkrs))
  
  p5 <- outlier_df %>%
    mutate(flag = ratio > 2) %>%
    ggplot(aes(x = reorder(marker, ratio),
               y = ratio, fill = flag)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("FALSE" = "#1D9E75","TRUE" = "#D85A30"),
                      name = "Ratio > 2") +
    geom_hline(yintercept = 2, linetype = "dashed", colour = "grey40") +
    coord_flip() +
    labs(title    = "Check 5: Per-marker 99th/95th percentile ratio",
         subtitle = "Orange bars = potential outlier-dominated markers",
         x = NULL, y = "p99 / p95 ratio") +
    theme_qc
  save_plot(p5, "05_outlier_ratios", w = 8, h = max(5, length(marker_cols) * 0.25 + 2))
  
  # ── 6. SAMPLE MIXING / BATCH CHECK ──────────────────────────────────────────
  cat("── Check 6/8: Sample mixing (PCA of sample medians)\n")
  
  med_mat <- expr %>%
    group_by(.data[[sid_col]]) %>%
    summarise(across(all_of(marker_cols), median, na.rm = TRUE),
              .groups = "drop")
  
  pca_res  <- prcomp(med_mat[, marker_cols], scale. = TRUE)
  pca_df   <- as.data.frame(pca_res$x[, 1:2])
  pca_df$sample    <- med_mat[[sid_col]]
  pca_df$condition <- metadata[pca_df$sample, "condition"]
  pct_var  <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  
  # Rough check: do samples from same condition cluster more tightly than across?
  if (length(unique(pca_df$condition)) >= 2 && nrow(pca_df) >= 4) {
    d_mat   <- dist(pca_df[, c("PC1","PC2")])
    hc      <- hclust(d_mat)
    cond_v  <- pca_df$condition[hc$order]
    runs    <- rle(cond_v)$lengths
    status_6 <- if (max(runs) == length(cond_v)) "warn" else "pass"
    flag("batch_mixing", status_6,
         if (status_6 == "warn")
           "Samples separate completely by condition in PCA — consider whether this is biology or batch"
         else
           "Samples mix reasonably across conditions in PCA space")
  } else {
    status_6 <- "pass"
    flag("batch_mixing", "pass", "Too few samples for batch separation test")
  }
  
  p6 <- ggplot(pca_df, aes(x = PC1, y = PC2,
                           colour = condition, label = sample)) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    labs(title    = "Check 6: PCA of per-sample median expression",
         subtitle = "Separation by sample > condition suggests batch effect",
         x = paste0("PC1 (", pct_var[1], "%)"),
         y = paste0("PC2 (", pct_var[2], "%)")) +
    theme_qc
  save_plot(p6, "06_pca_sample_mixing", w = 7, h = 5)
  
  # ── 7. MARKER CORRELATION (unmixing check) ───────────────────────────────────
  cat("── Check 7/8: Marker correlation structure\n")
  
  sub10k  <- expr %>% slice_sample(n = min(10000, nrow(expr)))
  cor_mat <- cor(sub10k[, marker_cols], method = "spearman", use = "pairwise.complete.obs")
  
  # Check mutual exclusion pairs if present
  me_pairs <- list(
    c("CD4","CD8"), c("CD3","CD19"), c("CD4","CD14"), c("CD3","CD14")
  )
  me_results <- lapply(me_pairs, function(pair) {
    cols <- marker_cols[grepl(paste0("^", pair[1], "|^", pair[2]), marker_cols)]
    if (length(cols) < 2) return(NULL)
    r <- cor(sub10k[[cols[1]]], sub10k[[cols[2]]],
             method = "spearman", use = "complete.obs")
    data.frame(pair = paste(pair, collapse = " vs "), r = round(r, 3),
               flag = r > 0.2)
  })
  me_df      <- do.call(rbind, Filter(Negate(is.null), me_results))
  n_bad_me   <- if (!is.null(me_df)) sum(me_df$flag) else 0
  status_7   <- if (n_bad_me > 0) "warn" else "pass"
  flag("marker_correlation", status_7,
       if (!is.null(me_df))
         paste("Mutual exclusion check:", paste(apply(me_df, 1, function(r)
           sprintf("%s r=%.2f%s", r["pair"], as.numeric(r["r"]),
                   if (as.logical(r["flag"])) " [!]" else "")), collapse="; "))
       else "No mutual exclusion pairs found to check")
  
  # Save correlation plot as PNG directly (corrplot doesn't return ggplot)
  cor_path <- file.path(plot_dir, "07_marker_correlation.png")
  png(cor_path, width = 1400, height = 1300, res = 150)
  corrplot::corrplot(cor_mat,
                     method = "color", type = "full", order = "hclust",
                     hclust.method = "ward.D2",
                     tl.cex = 0.6, tl.col = "black",
                     col = colorRampPalette(c("#185FA5","white","#A32D2D"))(200),
                     title = "Check 7: Marker correlation structure (Spearman)",
                     mar = c(0, 0, 2, 0))
  dev.off()
  plots[["07_marker_correlation"]] <- cor_path
  
  # ── 8. NAIVE CELL TYPE PROPORTIONS ──────────────────────────────────────────
  cat("── Check 8/8: Naive cell type proportions\n")
  
  # Detect lineage marker columns by partial name
  find_marker <- function(patterns, cols) {
    for (p in patterns) {
      m <- cols[grepl(p, cols, ignore.case = TRUE)]
      if (length(m) > 0) return(m[1])
    }
    return(NULL)
  }
  
  cd3_col  <- find_marker(c("CD3"),  marker_cols)
  cd4_col  <- find_marker(c("CD4"),  marker_cols)
  cd8_col  <- find_marker(c("CD8"),  marker_cols)
  cd19_col <- find_marker(c("CD19"), marker_cols)
  cd14_col <- find_marker(c("CD14"), marker_cols)
  cd56_col <- find_marker(c("CD56","NKp46"), marker_cols)
  
  found_lineage <- !sapply(list(cd3_col,cd4_col,cd8_col,cd19_col,cd14_col), is.null)
  cat(sprintf("  Lineage markers found: CD3=%s CD4=%s CD8=%s CD19=%s CD14=%s\n",
              cd3_col, cd4_col, cd8_col, cd19_col, cd14_col))
  
  if (sum(found_lineage) >= 3) {
    thr <- 1.5  # arcsinh positivity threshold
    
    prop_df <- expr %>%
      mutate(
        T_cell = if (!is.null(cd3_col))  .data[[cd3_col]]  > thr else FALSE,
        CD4T   = if (!is.null(cd4_col))  .data[[cd4_col]]  > thr & T_cell else FALSE,
        CD8T   = if (!is.null(cd8_col))  .data[[cd8_col]]  > thr & T_cell else FALSE,
        B_cell = if (!is.null(cd19_col)) .data[[cd19_col]] > thr else FALSE,
        Mono   = if (!is.null(cd14_col)) .data[[cd14_col]] > thr else FALSE,
        NK     = if (!is.null(cd56_col)) .data[[cd56_col]] > thr &
          (if (!is.null(cd3_col)) .data[[cd3_col]] < thr else TRUE) else FALSE
      ) %>%
      group_by(.data[[sid_col]], condition) %>%
      summarise(
        T_cell = mean(T_cell, na.rm=TRUE)*100,
        CD4T   = mean(CD4T,   na.rm=TRUE)*100,
        CD8T   = mean(CD8T,   na.rm=TRUE)*100,
        B_cell = mean(B_cell, na.rm=TRUE)*100,
        Mono   = mean(Mono,   na.rm=TRUE)*100,
        NK     = mean(NK,     na.rm=TRUE)*100,
        .groups = "drop"
      )
    
    # Flag samples where T cells < 20% or B cells > 30%
    n_bad_prop <- sum(prop_df$T_cell < 20 | prop_df$B_cell > 30, na.rm=TRUE)
    status_8   <- if (n_bad_prop > 0) "warn" else "pass"
    flag("cell_proportions", status_8,
         sprintf("%d samples with unusual lineage proportions (T<20%% or B>30%%)",
                 n_bad_prop))
    
    long_prop <- prop_df %>%
      pivot_longer(c(T_cell,CD4T,CD8T,B_cell,Mono,NK),
                   names_to="population", values_to="pct")
    
    p8 <- ggplot(long_prop, aes(x = .data[[sid_col]],
                                y = pct, fill = population)) +
      geom_bar(stat = "identity") +
      facet_wrap(~condition, scales = "free_x") +
      coord_flip() +
      labs(title    = "Check 8: Naive cell type proportions (arcsinh threshold = 1.5)",
           subtitle = "Rough estimate only — not a replacement for proper gating",
           x = NULL, y = "% of cells") +
      scale_fill_manual(values = c(
        T_cell = "#185FA5", CD4T = "#5DCAA5", CD8T = "#1D9E75",
        B_cell = "#7F77DD", Mono = "#D85A30",  NK   = "#EF9F27")) +
      theme_qc
    save_plot(p8, "08_cell_proportions", w = 8, h = 5)
  } else {
    flag("cell_proportions","warn","Too few lineage markers identified for proportion check")
    prop_df <- NULL
  }
  
  # ── GENERATE HTML REPORT ─────────────────────────────────────────────────────
  cat("\n── Generating HTML report\n")
  
  status_colours <- c(pass="#1D9E75", warn="#BA7517", fail="#A32D2D")
  status_bg      <- c(pass="#E1F5EE", warn="#FAEEDA", fail="#FCEBEB")
  
  summary_rows <- paste(sapply(names(msgs), function(nm) {
    m   <- msgs[[nm]]
    col <- status_colours[m$status]
    bg  <- status_bg[m$status]
    sprintf(
      '<tr style="background:%s">
         <td style="padding:8px 12px;font-weight:500;color:%s">%s</td>
         <td style="padding:8px 12px;color:%s">%s</td>
         <td style="padding:8px 12px;color:#444">%s</td>
       </tr>',
      bg, col, m$icon,
      col, gsub("_"," ", nm),
      m$msg)
  }), collapse="\n")
  
  n_pass <- sum(sapply(msgs, function(m) m$status=="pass"))
  n_warn <- sum(sapply(msgs, function(m) m$status=="warn"))
  n_fail <- sum(sapply(msgs, function(m) m$status=="fail"))
  
  plot_sections <- paste(sapply(names(plots), function(nm) {
    rel_path <- file.path("plots", basename(plots[[nm]]))
    label    <- gsub("_"," ", gsub("^[0-9]+_","", nm))
    sprintf('<div style="margin:2rem 0">
               <h3 style="font-size:15px;font-weight:500;margin-bottom:8px">%s</h3>
               <img src="%s" style="max-width:100%%;border-radius:8px;
                    border:1px solid #e0e0e0" alt="%s">
             </div>', label, rel_path, label)
  }), collapse="\n")
  
  html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>cyCONDOR QC Report</title>
<style>
  body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
       max-width:900px;margin:2rem auto;padding:0 1.5rem;
       color:#1a1a1a;background:#fafafa;line-height:1.6}
  h1{font-size:22px;font-weight:600;margin-bottom:4px}
  h2{font-size:17px;font-weight:500;margin:2.5rem 0 1rem;
     border-bottom:1px solid #e0e0e0;padding-bottom:6px}
  .meta{font-size:13px;color:#666;margin-bottom:2rem}
  .summary-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin:1.5rem 0}
  .stat-card{padding:1rem;border-radius:8px;text-align:center}
  .stat-card .n{font-size:28px;font-weight:600}
  .stat-card .l{font-size:12px;margin-top:2px}
  .pass-card{background:#E1F5EE;color:#0F6E56}
  .warn-card{background:#FAEEDA;color:#633806}
  .fail-card{background:#FCEBEB;color:#791F1F}
  table{width:100%;border-collapse:collapse;font-size:14px;margin-top:8px}
  th{text-align:left;padding:8px 12px;font-size:12px;font-weight:600;
     text-transform:uppercase;letter-spacing:.05em;
     background:#f0f0f0;border-bottom:2px solid #ddd}
  tr+tr td{border-top:1px solid #eee}
  .footer{font-size:12px;color:#999;margin-top:3rem;border-top:1px solid #eee;padding-top:1rem}
</style>
</head>
<body>

<h1>cyCONDOR post-load QC report</h1>
<div class="meta">
  Generated: %s &nbsp;|&nbsp;
  Cells: %s &nbsp;|&nbsp;
  Markers: %d &nbsp;|&nbsp;
  Samples: %d &nbsp;|&nbsp;
  arcsinh cofactor: %g
</div>

<div class="summary-grid">
  <div class="stat-card pass-card"><div class="n">%d</div><div class="l">Passed</div></div>
  <div class="stat-card warn-card"><div class="n">%d</div><div class="l">Warnings</div></div>
  <div class="stat-card fail-card"><div class="n">%d</div><div class="l">Failed</div></div>
</div>

<h2>Check summary</h2>
<table>
  <tr><th>Status</th><th>Check</th><th>Details</th></tr>
  %s
</table>

<h2>Plots</h2>
%s

<div class="footer">
  Run with cyCONDOR QC wrapper v1.0 &nbsp;|&nbsp;
  Adjust marker name patterns in find_marker() if lineage detection failed.
</div>

</body>
</html>',
                  format(Sys.time(), "%Y-%m-%d %H:%M"),
                  format(nrow(expr), big.mark=","),
                  length(marker_cols),
                  length(unique(expr[[sid_col]])),
                  cofactor,
                  n_pass, n_warn, n_fail,
                  summary_rows,
                  plot_sections
  )
  
  report_path <- file.path(out_dir, "condor_qc_report.html")
  writeLines(html, report_path)
  
  cat(sprintf("\n══ QC complete: %d passed, %d warnings, %d failed\n",
              n_pass, n_warn, n_fail))
  cat(sprintf("   Report saved to: %s\n\n", report_path))
  
  invisible(list(
    summary  = msgs,
    plots    = plots,
    report   = report_path,
    expr_qc  = expr,        # expression matrix with metadata joined
    outliers = outlier_df,
    cors     = cor_mat
  ))
}

# null-coalescing helper (base R doesn't have %||%)
`%||%` <- function(a, b) if (!is.null(a)) a else b
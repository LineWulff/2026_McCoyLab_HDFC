
plot_pseudobulkPCA <- function(condor, extra_anno = NULL, marker_cols = colnames(condor$expr$orig), AnnoTable=AnnoTable){
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

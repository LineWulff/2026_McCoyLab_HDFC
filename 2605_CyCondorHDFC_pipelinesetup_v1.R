#' R script CyCondor pipeline for high dimensional flow cytometry data
#' Author: Line Wulff
#' Using: https://lorenzobonaguro.github.io/cyCONDOR/articles/cyCONDOR.html
#' Date (created): 26-05-12
rm(list=ls())

# working directory, where you want the analysis to reside, should include all R scripts 
# get folder path for samples and working directory and copy in
# mac: right click folder and press option, copy "" as pathname
# windows: when in folder click top panel with folder path and copy
# make sure path ends with /
main_dir <- "/Users/linewulff/Documents/work/projects/2026_McCoyLab_HDFC/"
setwd(main_dir)

#### ---- libraries used - MANDATORY ---- ####
library(cyCONDOR)
library(tidyr)
library(tidyverse)
library(stringr)
library(ggrastr)
library(ggplot2)
library(reshape)
library(pheatmap)
library(uwot)
library(viridis)
library(scales)
source("find_high_cor_pairs.R") # should be saved in working directory
source("plot_pseudobulkPCA.R") # -||-
source("perc_dist.R") # -||-
source("geom_cluster_labels.R")  # -||-

#### ---- Variables to use for script - MANDATORY--- ####
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
proj <- "CovidHealthPBMC" #project name, will be on all plots

#### ---- Read in data - MANDATORY ---- ####
# get folder path for samples and working directory and copy in
# mac: right click folder and press option, copy "" as pathname
# windows: when in folder click top panel with folder path and copy
# make sure path ends with /
sample_dir <- "/Users/linewulff/Documents/work/projects/2026_McCoyLab_HDFC/test_data/27940719/"
# in folder need csv file named AnnoTable, check example in
# /Promise RAID/Journal Club/Line/HDC_wCyCondor

condor <- prep_fcd(data_path = sample_dir, 
                   cross_path_with_anno = TRUE,
                   max_cell = 10000000, # 10000000 start w. high limit to include all
                   useCSV = TRUE, # change to false if using fcs files
                   transformation = "auto_logi", #"arcsinh",
                   #remove_param = c("FSC-H", "SSC-H", "FSC-W", "SSC-W", "Time"), 
                   anno_table = paste0(sample_dir,"AnnoTable.csv"), 
                   filename_col = "filename")

AnnoTable <- read.csv(paste0(sample_dir,"AnnoTable.csv"), row.names = 1)

# all meta data stored per cell in below, you can add, change and color/split plots/stats by anything added here as you go
head(condor$anno$cell_anno)
# Check that your AnnoTable per sample looks correct
AnnoTable


#### ---- Quality controls - not included in CyCondor pipeline - MANDATORY ---- ####
# makes a QC folder where plots of below QC measures and 
if (file.exists("QC")){} else {dir.create(file.path(main_dir, "QC"))}

# check how many cells are attributed from each sample and go back to prep_fcd and change max_cell
# to minimum of below to avoid downstream effects of larger sample sizes vs smaller sample sizes
# condition groups should contribute proprotinally
cell_no <- as.data.frame(table(condor$anno$cell_anno$sample_id))
colnames(cell_no) <- c("sample_id","cell_no")
cell_no

ggplot(cell_no, aes(x=sample_id, y=cell_no)) +
  geom_bar(stat="identity")+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90))+
  labs(y="cell number per sample")+
  geom_hline(yintercept = min(cell_no$cell_no), linetype="dashed")
ggsave(paste0(main_dir,"QC/",dato,proj,"_CellNumberPerSample.pdf"), width = 5, height = 5, units = "in")

# percentage distribution between groups
perc_dist(condor, anno = "condition")
perc_dist(condor, anno = "sex")


# data frame for QC purposes
df_long <- melt(cbind(as.data.frame(condor$expr$orig),
                      sample_id=condor$anno$cell_anno$sample_id,
                      condition=condor$anno$cell_anno$condition), 
                id.vars = c("sample_id","condition"))

# flourescence as histogram per marker - transformation QC
# Check bimodal / unimodal (markers for rare cell type markers) / trimodal (markers w. neg, low and high exp patterns)
ggplot(df_long, aes(value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()
ggsave(paste0(main_dir,"QC/",dato,proj,"_QC-FlourescenceHistograms.pdf"), width = 10, height = 10, units = "in")

# Flourescence intensity per marker per sample - as in original FlowSOM paper
# check samples look similar per condition
ggplot(df_long, aes(x=sample_id, y=value, color = condition))+
  geom_boxplot()+
  theme_classic()+
  facet_wrap(~variable, scales = "free")+
  theme(axis.text.x = element_text(angle=90))
ggsave(paste0(main_dir,"QC/",dato,proj,"_QC-FlourescenceBoxPlotsSplitPerSample.pdf"), width = 15, height = 15, units = "in")


# 99% saturation check
# Should be 0.001, if any higher inspect marker individually
apply(as.data.frame(condor$expr$orig), 2, function(x) mean(x > quantile(x, 0.999)))

# Correlation of markers, quick check for spillovers/potential compensation issues
# check spill over/real correlation of markers in heatmap
cor_mat <- cor(condor$expr$orig)
pheatmap(cor_mat)
ggsave(paste0(main_dir,"QC/",dato,proj,"_QC-FlourescenceCorrelation.pdf"), width = 8, height = 8, units = "in")

# function to identify high correlation pairs
cor_pairs <- find_high_cor_pairs(cor_mat, threshold = 0.2) 
cor_pairs # all pairs with cor > 0.2 should be checked 
# If biologically sound can be included, e.g. T cell markers like CD4 and CD3 can correlate,
# Mutually exclusive markers like Cd4 and Cd8 should not correlate
# If there's significant spill over you cannot use those markers for dimensionality reduction and clustering
# Algorithm cannot know this is not true signal 

# PCA on per-sample medians before UMAP (check for batch effects/normal variance between samples)
marker_cols <- colnames(condor$expr$orig) # or choose specific ones
# extra_anno accepts up to two arguments, if you don't want to color or shape by anything set extra_anno = NULL
PCA1 <- plot_pseudobulkPCA(condor, extra_anno = c("condition","sex"), marker_cols = marker_cols, AnnoTable=AnnoTable)
print(PCA1)
ggsave(paste0(main_dir,"QC/",dato,proj,"_QC-PseudobulkPCA.pdf"), width = 6, height = 6, units = "in")

#### ---- Batch correction ---- ####
# If you run UMAP and clustering and clusters do not make biologically sense 
# They could be affected by batch effects
# Run this section and retry
condor <- runPCA(fcd = condor, 
                 data_slot = "orig")

condor <- harmonize_PCA(fcd = condor, 
                        batch_var = c("sample_id"),
                        data_slot = "orig")

condor <- runUMAP(fcd = condor, 
                  input_type = "pca", 
                  data_slot = "norm",
                  prefix= NULL)


#### ---- Dimensionality reduction - WITHOUT batch correction ---- ####
# condor standard runs umap/tsne on PCs, however, not advisable unless you have strong batch effects.
# check that marker cols contains the markers you want to calculate umap and clustering on
marker_cols
# run umap with uwot package - takes some time!!
umap_emb <- umap(condor$expr$orig[,marker_cols])
condor$umap$orig <- umap_emb[,c(1,2)]; colnames(condor$umap$orig) <- c("UMAP1","UMAP2") 

#### ---- Plotting UMAPs ---####
# first set whether you're using batch corrected or not
batchcor <- "Yes" # "Yes" or "No"
## first creating a folder
if(batchcor=="No"){if (file.exists("UMAPs")){} else {dir.create(file.path(main_dir, "UMAPs"))};umap <-"orig";umap_dir <- "UMAPs/";input_expr <- c("expr","orig")} else if(batchcor=="Yes"){if (file.exists("UMAPs_bc")){} else {dir.create(file.path(main_dir, "UMAPs_bc"))};umap <-"pca_norm";umap_dir <- "UMAPs_bc/";input_expr <- c("pca","norm")}

# color by condition
ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$condition))+
  geom_point_rast()+
  labs(color="condition")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAP_conditions.pdf"), width = 5, height = 5, units = "in")

# density version of above
ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$condition))+
  geom_density_2d(bins=20)+
  labs(color="condition")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAP_conditioncontour.pdf"), width = 5, height = 5, units = "in")

ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$condition))+
  geom_point_rast()+
  facet_wrap(~condor$anno$cell_anno$condition, scales = "free")+
  labs(color="condition")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAPsplit_condition.pdf"), width = 10, height = 5, units = "in")

#color by sampleid
ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sample_id))+
  geom_point_rast()+
  labs(color="sample_id")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAP_sampleid.pdf"), width = 5, height = 5, units = "in")

# density version of above
ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sample_id))+
  geom_density_2d(bins=20)+
  labs(color="sample_id")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAP_sampleidcontour.pdf"), width = 5, height = 5, units = "in")


# color by sex
ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sex))+
  geom_point_rast()+
  labs(color="sex")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAP_sex.pdf"), width = 5, height = 5, units = "in")

# split and colored by sex
ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sex))+
  geom_point_rast()+
  facet_wrap(~condor$anno$cell_anno$sex, scales = "free")+
  labs(color="sex")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAPsplit_sex.pdf"), width = 10, height = 5, units = "in")

# color by sex split by condition
ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sex))+
  geom_point_rast()+facet_wrap(~condor$anno$cell_anno$condition, scales = "free")+
  labs(color="sex")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAP_conditionXsex.pdf"), width = 10, height = 5, units = "in")


# plot UMAP of each marker expression
for (mark in marker_cols){
  pdf(paste0(main_dir,umap_dir,dato,proj,"_UMAP_",mark,".pdf"), width = 5, height = 5)
  UMAP1 <- ggplot(condor$umap[[umap]], 
         aes(x=UMAP1, y=UMAP2, color = condor$expr$orig[,mark]))+
    geom_point_rast()+
    scale_colour_viridis_c()+
    labs(color=mark)+
    theme_classic()+
    theme(axis.ticks = element_blank(),axis.text = element_blank())
  print(UMAP1)
  dev.off()
}


#### ---- Clustering w. FlowSOM ---- ####
# first set whether you're using batch corrected or not
batchcor <- "Yes" # "Yes" or "No"
## first creating a folder
if(batchcor=="No"){if (file.exists("UMAPs")){} else {dir.create(file.path(main_dir, "UMAPs"))};umap <-"orig";umap_dir <- "UMAPs/";input_expr <- c("expr","orig")} else if(batchcor=="Yes"){if (file.exists("UMAPs_bc")){} else {dir.create(file.path(main_dir, "UMAPs_bc"))};umap <-"pca_norm";umap_dir <- "UMAPs_bc/";input_expr <- c("pca","norm")}


# run multiple resolutions/cluster numbers for flowsom at once
clus_res <- seq(5,18) # calculating from 5 to 15 clusters
for (res in clus_res){
  print(res)
  # calculating
  condor <- runFlowSOM(fcd = condor, 
                     input_type = input_expr[1], 
                     data_slot = input_expr[2], 
                     nClusters = res)}

# or run just one by removing all in front of code #
# condor <- runFlowSOM(fcd = condor, 
#                      input_type = input_expr, 
#                      data_slot = "orig", 
#                      nClusters = 10)
clus_res <- seq(20,60,5)
for (res in clus_res){
  print(res)
  # calculating
  condor <- runPhenograph(fcd = condor, 
                          input_type = input_expr[1], 
                          data_slot = input_expr[2], 
                          k = res)}

#### ---- UMAPs of all available clusterings ---- ####
# run loop to save all your clusterings in their respective umap folders
for (clus in names(condor$clustering)){
  # is it batch corrected?
  if (str_detect(clus,"norm")){umap_dir <- "UMAPs_bc/";umap <-"pca_norm"} else {umap_dir <- "UMAPs/";umap <-"orig"}
  # Detect different clustering methods
  if (startsWith(clus,"Flow")){clus_var = "FlowSOM"} else {clus_var = "phenograph"}
  if (length(unique(condor$clustering[[clus]][,1]))<50){
    plot_df <- cbind(as.data.frame(condor$umap[[umap]]),clustering=condor$clustering[[clus]][,1])
    
  p1 <- ggplot(plot_df, 
               aes(x=UMAP1, y=UMAP2, color = clustering))+
    geom_point_rast()+
    labs(color="", title = clus)+
    theme_classic()+
    theme(axis.ticks = element_blank(),axis.text = element_blank())
  p1 <- p1+geom_cluster_labels(plot_df, "clustering")
  pdf(paste0(main_dir,umap_dir,dato,proj,"_UMAP_",clus,".pdf"), width = 7, height = 5)
  print(p1)
 dev.off()}
  else {print(paste(clus, "had 50 < clusters and was not plotted."))}
}

## now make sure your standard UMAP is set to your batch/not batch corrected again
# first set whether you're using batch corrected or not
batchcor <- "Yes" #or "No"
## first creating a folder
if(batchcor=="No"){if (file.exists("UMAPs")){} else {dir.create(file.path(main_dir, "UMAPs"))};umap <-"orig";umap_dir <- "UMAPs/";input_expr <- c("expr","orig")} else if(batchcor=="Yes"){if (file.exists("UMAPs_bc")){} else {dir.create(file.path(main_dir, "UMAPs_bc"))};umap <-"pca_norm";umap_dir <- "UMAPs_bc/";input_expr <- c("pca","norm")}

#### ---- Inspection of specific clusterings ---- ####
## color UMAP by cluster
# Set res to any of the clustering run below and plot for this resolution
# check avalable clusterings, res should match one of these
names(condor$clustering)
# check number of clusters in a particular phenograph clustering
length(unique(condor$clustering[["FlowSOM_pca_norm_k_14" ]][,1]))
## FlowSOM clusterings examples
res <- "FlowSOM_expr_orig_k_7" # NO batch correction
res <- "FlowSOM_pca_norm_k_14" # batch corrected
#Phenograph clustering example
res <- "phenograph_expr_orig_k_80" # NO natch correction
res <- "phenograph_pca_norm_k_34" # batch corrected
# run below line for ease of plotting further down
if (startsWith(res,"Flow")){clus_var = "FlowSOM"} else {clus_var = "Phenograph"}
  
# split per condition
p1 <- ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$clustering[[res]][,1]))+
  geom_point_rast()+
  facet_wrap(~condor$anno$cell_anno$condition, scales = "free")+
  labs(color=res)+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
print(p1)
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAPsplit_conditionx",res,".pdf"), width = 12, height = 5, units = "in", plot =p1)

# one UMAP
p1 <- ggplot(condor$umap[[umap]], 
       aes(x=UMAP1, y=UMAP2, color = condor$clustering[[res]][,1]))+
  geom_point_rast()+
  labs(color=res)+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
print(p1)
ggsave(paste0(main_dir,umap_dir,dato,proj,"_UMAP_",res,".pdf"), width = 7, height = 5, units = "in", plot =p1)

#### ---- Marker expression - RidgePlots per cluster and heatmap ---- ####
# first new folder per clustering
resk_col <- hue_pal()(length(unique(condor$clustering[[res]][,1]))); names(resk_col) <- seq(1,length(unique(condor$clustering[[res]][,1]))) 
if (file.exists(res)){} else {dir.create(file.path(main_dir, res))}

# Histogram of expression for each marker across different clusters
for (mark in marker_cols){
    histo_exp_plot <- ggplot(cbind(condor$expr$orig,cluster=condor$clustering[[res]][,1]),
                             aes(x=condor$expr$orig[,mark], fill=cluster))+
      geom_density()+
      facet_wrap(.~cluster, ncol = 1)+
      labs(fill=res, x=mark)+
      theme_classic()
    if (length(unique(condor$clustering[[res]][,1]))>5){len=1.1*length(unique(condor$clustering[[res]][,1]))}else(len=8)
    ggsave(paste0(main_dir,res,"/",dato,proj,"_Histogram_",mark,".pdf"), width = 5, height = len, plot = histo_exp_plot)}

# Histogram of expression for each marker across different clusters, split between conditions
# change "condition" to any (non continous) column name from your meta data - change in the two marked spots
for (mark in marker_cols){
  histo_exp_plot <- ggplot(cbind(condor$expr$orig,cluster=condor$clustering[[res]][,1],condor$anno$cell_anno),
                           aes(x=condor$expr$orig[,mark], fill = cluster))+
    geom_density()+
    facet_grid(cluster~condition)+ # change ~????
    labs(fill=res, x=mark)+
    theme_classic()
  if (length(unique(condor$clustering[[res]][,1]))>5){len=1.1*length(unique(condor$clustering[[res]][,1]))}else(len=8)
  ggsave(paste0(main_dir,res,"/",dato,proj,"_Histogram_",mark,"conditionsplit.pdf"), width = 7, height = len, plot = histo_exp_plot)}

# Or profiles of the conditions overlaid each other
for (mark in marker_cols){
  histo_exp_plot <- ggplot(cbind(condor$expr$orig,cluster=condor$clustering[[res]][,1],condor$anno$cell_anno),
                           aes(x=condor$expr$orig[,mark], fill = condition, color=condition))+
    geom_density(alpha=0.1)+
    facet_wrap(.~cluster, ncol = 1)+
    labs(fill=res, x=mark)+
    guides(fill = "none")+
    theme_classic()
  if (length(unique(condor$clustering[[res]][,1]))>5){len=1.1*length(unique(condor$clustering[[res]][,1]))}else(len=8)
  ggsave(paste0(main_dir,res,"/",dato,proj,"_Histogram_",mark,"conditioncolor.pdf"), width = 4, height = len, plot = histo_exp_plot)}

# Heatmap of protein expression averaged per cluster
p1 <- plot_marker_HM(fcd = condor,
              cluster_rows = TRUE, cluster_cols = TRUE,
               expr_slot = "orig",
               cluster_slot = res,
               cluster_var = clus_var) 
print(p1)
ggsave(paste0(main_dir,res,"/",dato,proj,"_heatmap_avScaExp.pdf"), width = 6, height = 10, units = "in", plot = p1)

#### ---- Sample distribution plots - NOT DONE STOP HERE ---- ####
p1 <- plot_confusion_HM(fcd = condor,
                  cluster_slot = res, 
                  cluster_var = clus_var,
                  group_var = "condition", 
                  size = 30)
print(p1)
ggsave(paste0(main_dir,res,"/",dato,proj,"_confusionMatrix.pdf"), width = 0.5*length(unique(condor$clustering[[res]][,1])), height = 2, units = "in", plot = p1)

# barplot - dots per sample dist
p1 <- condor_cluster_composition(
  condor      = condor,
  cluster_slot = res,
  anno_col    = "condition",
  palette = c("Healthy"="#00BFC4","COVID19"="#F8766D"), # use default, or specify in this line, # in front for default
  metric      = "percentage", # count or percentage
  errorbar_fun = "sd")
print(p1)
ggsave(paste0(main_dir,res,"/",dato,proj,"_SampleDistribution.pdf"), width = 0.33*length(unique(condor$clustering[[res]][,1])), height = 5, units = "in", plot = p1)

#### ---- Label and subset condor object - metaclustering ---- #### 
# replace and extend labels depending on your clustering
# you need labels for all cell, but can be NA/unkown
# here you can also merge vy having the same ID for multiple cluster numbers
condor <- metaclustering(fcd = condor, 
                         cluster_slot = res, 
                         cluster_var = clus_var, 
                         cluster_var_new = "metaclusters", 
                         metaclusters = c("1" = "Classical Monocytes",
                                          "2" = "CD4 CD45RA+ CD127+",
                                          "3" = "CD8 CD45RA+ CD127+", 
                                          "4" = "NK dim",
                                          "5" = "CD8 CD45RA+ CD127-",
                                          "6" = "Classical Monocytes",
                                          "7" = "Unconventional T cells", 
                                          "8" = "CD4 CD45RA- CD127+",
                                          "9" = "CD16+ Monocytes",
                                          "10" = "CD4 CD127-",
                                          "11" = "Classical Monocytes", 
                                          "12" = "CD8 CD45RA- CD127+", 
                                          "13" = "CD8 CD45RA- CD127+",
                                          "14" = "NK bright",
                                          "15" = "CD8 CD45RA+ CD127-",
                                          "16" = "CD4 CD25+",
                                          "17" = "B cells",
                                          "18" = "Unconventional T cells",
                                          "19" = "Classical Monocytes",
                                          "20" = "pDCs",
                                          "21" = "CD8 CD45RA+ CD127+",
                                          "22" = "Basophils",
                                          "23" = "Mixed",
                                          "24" = "B cells",
                                          "25" = "NK bright"))

## If you want to combine clusters from another clustering
condor$clustering[["metaclusters"]] <- condor$clustering[[res]]
# pick cells from the other clustering
rare_cells <- rownames(condor$clustering[["phenograph_pca_norm_k_25"]][condor$clustering[["phenograph_pca_norm_k_25"]][,"Phenograph"] == 21,])
condor$clustering[["metaclusters"]][rare_cells,"metaclusters"] <- "pDC"


# you can now set your res to your annotated clustering and rerun your plots from above
res <- "annotated_clustering"

#### ----####

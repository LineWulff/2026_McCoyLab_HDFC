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
                   useCSV = TRUE, 
                   transformation = "auto_logi", #"arcsinh",
                   #remove_param = c("FSC-H", "SSC-H", "FSC-W", "SSC-W", "Time"), 
                   anno_table = paste0(sample_dir,"AnnoTable_ext.csv"), 
                   filename_col = "filename")

AnnoTable <- read.csv(paste0(sample_dir,"AnnoTable_ext.csv"), row.names = 1)

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
pdf(paste0(main_dir,"QC/",dato,proj,"_CellNumberPerSample.pdf"), width = 5, height = 5)
ggplot(cell_no, aes(x=sample_id, y=cell_no)) +
  geom_bar(stat="identity")+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90))+
  labs(y="cell number per sample")+
  geom_hline(yintercept = min(cell_no$cell_no), linetype="dashed")
dev.off()
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
pdf(paste0(main_dir,"QC/",dato,proj,"_QC-FlourescenceHistograms.pdf"), width = 10, height = 10)
ggplot(df_long, aes(value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()
dev.off()

# Flourescence intensity per marker per sample - as in original FlowSOM paper
# check samples look similar per condition
pdf(paste0(main_dir,"QC/",dato,proj,"_QC-FlourescenceBoxPlotsSplitPerSample.pdf"), width = 15, height = 15)
ggplot(df_long, aes(x=sample_id, y=value, color = condition))+
  geom_boxplot()+
  theme_classic()+
  facet_wrap(~variable, scales = "free")+
  theme(axis.text.x = element_text(angle=90))
dev.off()

# 99% saturation check
# Should be 0.001, if any higher inspect marker individually
apply(as.data.frame(condor$expr$orig), 2, function(x) mean(x > quantile(x, 0.999)))

# Correlation of markers, quick check for spillovers/potential compensation issues
# check spill over/real correlation of markers in heatmap
cor_mat <- cor(condor$expr$orig)
pdf(paste0(main_dir,"QC/",dato,proj,"_QC-FlourescenceCorrelation.pdf"), width = 8, height = 8)
pheatmap(cor_mat)
dev.off()
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
pdf(paste0(main_dir,"QC/",dato,proj,"_QC-PseudobulkPCA.pdf"), width = 6, height = 6)
print(PCA1)
dev.off()

### ---- Dimensionality reduction ---- ####
if (file.exists("UMAPs")){} else {dir.create(file.path(main_dir, "UMAPs"))}
# condor standard runs umap/tsne on PCs, however, not advisable.
# check that marker cols contains the markers you want to calcualte umap and clustering on
marker_cols
# run umap with uwot package - takes some time!!
umap_emb <- umap(condor$expr$orig[,marker_cols])
condor$umap$orig <- umap_emb[,c(1,2)]; colnames(condor$umap$orig) <- c("UMAP1","UMAP2") 

# color by condition
pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAP_conditions.pdf"), width = 5, height = 5)
ggplot(condor$umap$orig, 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$condition))+
  geom_point_rast()+
  labs(color="condition")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAPsplit_condition.pdf"), width = 10, height = 5)
ggplot(condor$umap$orig, 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$condition))+
  geom_point_rast()+
  facet_wrap(~condor$anno$cell_anno$condition, scales = "free")+
  labs(color="condition")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

#color by sampleid
pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAP_sampleid.pdf"), width = 5, height = 5)
ggplot(condor$umap$orig, 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sample_id))+
  geom_point_rast()+
  labs(color="sample_id")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

# color by sex
pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAP_sex.pdf"), width = 5, height = 5)
ggplot(condor$umap$orig, 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sex))+
  geom_point_rast()+
  labs(color="sex")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAPsplit_sex.pdf"), width = 10, height = 5)
ggplot(condor$umap$orig, 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sex))+
  geom_point_rast()+
  facet_wrap(~condor$anno$cell_anno$sex, scales = "free")+
  labs(color="sex")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

# color by sex split by condition
pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAP_conditionXsex.pdf"), width = 10, height = 5)
ggplot(condor$umap$orig, 
       aes(x=UMAP1, y=UMAP2, color = condor$anno$cell_anno$sex))+
  geom_point_rast()+facet_wrap(~condor$anno$cell_anno$condition, scales = "free")+
  labs(color="sex")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

# plot UMAP of each marker expression
for (mark in marker_cols){
  pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAP_",mark,".pdf"), width = 5, height = 5)
  UMAP1 <- ggplot(condor$umap$orig, 
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
# run multiple resolutions/cluster numbers for flwosom at once
clus_res <- seq(5,15) # calculating from 5 to 15 clusters
for (res in clus_res){
  print(res)
  # calculating
  condor <- runFlowSOM(fcd = condor, 
                     input_type = "expr", 
                     data_slot = "orig", 
                     nClusters = res)}

# or run just one by removing all in front of code #
# condor <- runFlowSOM(fcd = condor, 
#                      input_type = "expr", 
#                      data_slot = "orig", 
#                      nClusters = 10)

## color UMAP by cluster
# Set res to any of the clustering run below and plot for this resolution
res <- 10

# split per condition
pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAPsplit_conditionxFlowSOMk",res,".pdf"), width = 10, height = 5)
ggplot(condor$umap$orig, 
       aes(x=UMAP1, y=UMAP2, color = condor$clustering[[res]]$FlowSOM))+
  geom_point_rast()+
  facet_wrap(~condor$anno$cell_anno$condition, scales = "free")+
  labs(color=paste0("FlowSOM_k", res))+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

# one UMAP
pdf(paste0(main_dir,"UMAPs/",dato,proj,"_UMAPsplit_FlowSOMk",res,".pdf"), width = 5, height = 5)
ggplot(condor$umap$orig, 
       aes(x=UMAP1, y=UMAP2, color = condor$clustering[[res]]$FlowSOM))+
  geom_point_rast()+
  labs(color=paste0("FlowSOM_k", res))+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

#### ---- Marker expression - RidgePlots per cluster and heatmap ---- ####
# first new folder per clustering
resk_col <- hue_pal()(res); names(resk_col) <- seq(1,res)
if (file.exists(paste0("FlowSOM_k",res))){} else {dir.create(file.path(main_dir, paste0("FlowSOM_k",res)))}

# Histogram of expression for each marker across different clusters
for (mark in marker_cols){
    pdf(paste0(main_dir,"FlowSOM_k",res,"/",dato,proj,"_UMAP_",mark,".pdf"), width = 4, height = 8)
    histo_exp_plot <- ggplot(cbind(condor$expr$orig,FlowSOM=condor$clustering[[res]]$FlowSOM),
                             aes(x=condor$expr$orig[,mark], fill=FlowSOM))+
      geom_density()+
      facet_wrap(.~FlowSOM, ncol = 1)+
      labs(fill=paste0("FlowSOM k",res), x=mark)+
      theme_classic()
    print(histo_exp_plot)
    dev.off()}

# Histogram of expression for each marker across different clusters, split between conditions
# change "condition" to any (non continous) column name from your meta data - change in the two marked spots
for (mark in marker_cols){
  pdf(paste0(main_dir,"FlowSOM_k",res,"/",dato,proj,"_UMAP_",mark,"conditioncolor.pdf"), width = 4, height = 8)
  histo_exp_plot <- ggplot(cbind(condor$expr$orig,FlowSOM=condor$clustering[[res]]$FlowSOM,condor$anno$cell_anno),
                           aes(x=condor$expr$orig[,mark], fill = FlowSOM))+
    geom_density()+
    facet_grid(FlowSOM~condition)+ # change ~????
    labs(fill=paste0("FlowSOM k",res), x=mark)+
    theme_classic()
  print(histo_exp_plot)
  dev.off()}

# Or profiles of the conditions overlaid each other
for (mark in marker_cols){
  pdf(paste0(main_dir,"FlowSOM_k",res,"/",dato,proj,"_UMAP_",mark,"conditioncolor.pdf"), width = 4, height = 8)
  histo_exp_plot <- ggplot(cbind(condor$expr$orig,FlowSOM=condor$clustering[[res]]$FlowSOM,condor$anno$cell_anno),
                           aes(x=condor$expr$orig[,mark], fill = condition, color=condition))+
    geom_density(alpha=0.1)+
    facet_wrap(.~FlowSOM, ncol = 1)+
    labs(fill=paste0("FlowSOM k",res), x=mark)+
    guides(fill = "none")+
    theme_classic()
  print(histo_exp_plot)
  dev.off()}

# Heatmap of protein expression averaged per cluster
plot_marker_HM(fcd = condor,
               expr_slot = "orig",
               cluster_slot = paste0("FlowSOM_expr_orig_k_",res),
               cluster_var = "FlowSOM")

#### ---- Sample distribution plots ---- ####
plot_confusion_HM(fcd = condor,
                  cluster_slot = paste0("FlowSOM_expr_orig_k_",res), 
                  cluster_var = "FlowSOM",
                  group_var = "condition", 
                  size = 30)

# barplot - dots per sample dist


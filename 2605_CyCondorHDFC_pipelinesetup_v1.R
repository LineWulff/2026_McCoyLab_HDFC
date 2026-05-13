#' R script CyCondor pipeline for high dimensional flow cytometry data
#' Author: Line Wulff
#' Using: https://lorenzobonaguro.github.io/cyCONDOR/articles/cyCONDOR.html
#' Date (created): 26-05-12
rm(list=ls())

#### ---- libraries used - MANDATORY ---- ####
library(cyCONDOR)
library(tidyr)
library(tidyverse)
library(ggrastr)
library(ggplot2)
library(reshape)
library(pheatmap)
source('QC_CyCondor.R')
source("find_high_cor_pairs.R")

#### ---- Variables to use for script - MANDATORY--- ####
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
proj <- "CovidHealthPBMC" #project name, will be on all plots

#### ---- Read in data - MANDATORY ---- ####
# get folder path for samples and working directory and copy in
# mac: right click folder and press option, copy "" as pathname
# windows: when in folder click top panel with folder path and copy
# make sure path ends with /
sample_dir <- "/Users/linewulff/Documents/work/projects/2026_McCoyLab_HDFC/test_data/27940719/"

# working directory, where you want the analysis to reside
main_dir <- "/Users/linewulff/Documents/work/projects/2026_McCoyLab_HDFC/"
setwd(main_dir)

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
# untransformed data will be deleted, but necessary for QC
condor_unt <- prep_fcd(data_path = sample_dir, 
                       cross_path_with_anno = TRUE,
                       max_cell = 10000000, # 10000000 start w. high limit to include all
                       useCSV = TRUE, 
                       transformation = "none",
                       #remove_param = c("FSC-H", "SSC-H", "FSC-W", "SSC-W", "Time"), 
                       anno_table = paste0(sample_dir,"AnnoTable_ext.csv"), 
                       filename_col = "filename")
AnnoTable <- read.csv(paste0(sample_dir,"AnnoTable_ext.csv"), row.names = 1)

# all meta data stored per cell in below, you can add, change and color/split plots/stats by anything added here as you go
head(condor$anno$cell_anno)

# check how many cells are attributed from each sample and go back to prep_fcd and change max_cell
# to minimum of below to avoid downstream effects of larger sample sizes vs smaller sample sizes
cell_no <- as.data.frame(table(condor$anno$cell_anno$sample_id))
colnames(cell_no) <- c("sample_id","cell_no")
cell_no

#### ---- Quality controls - not included in CyCondor pipeline - MANDATORY ---- ####
# makes a QC folder where plots of below QC measures and 
if (file.exists("QC")){} else {dir.create(file.path(main_dir, "QC"))}

# check how many cells are attributed from each sample and go back to prep_fcd and change max_cell
# to minimum of below to avoid downstream effects of larger sample sizes vs smaller sample sizes
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
apply(as.data.frame(condor$expr), 2, function(x) mean(x > quantile(x, 0.999)))

# Correlation of markers, quick check for spillovers/potential compensation issues
# check spill over/real correlation of markers in heatmap
cor_mat <- cor(condor$expr$orig)
pdf(paste0(main_dir,"QC/",dato,proj,"_QC-FlourescenceCorrelation.pdf"), width = 8, height = 8)
pheatmap(cor_mat)
dev.off()
# function to identify high correlation pairs
find_high_cor_pairs(cor_mat, threshold = 0.5) 
# If there's significant spill over you cannot use those markers for dimensionality reduction and clustering
# Algorithm cannot know this is not true signal but just spill over

#### ---- Dimensionality reduction ---- ####
#### ---- ---- ####
#### ---- ---- ####

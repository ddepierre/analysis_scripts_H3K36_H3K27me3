######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# SCRIPT PLOT HEATMAP FIGURE 1 D E HEATMAP

######################################################################################################################################################
######################################################################################################################################################
# R version 3.4.2 (2017-09-28) -- "Short Summer"
# Copyright (C) 2017 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

######################################################################################################################################################
### LIBRARY
######################################################################################################################################################
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library("GenomicRanges")
library("ggplot2")
library("ggpubr")
library("gplots")
library("gsubfn")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("rtracklayer")
'%ni%' = Negate('%in%')

## define workdir
workdir = ""

######################################################################################################################################################
### LOAD FUNCTION
######################################################################################################################################################


# Load FUNCTION
source(paste0(workdir, "FIGURES/R_PLOT_FUNCTION/Script_HEATMAP_profile.R"))


######################################################################################################################################################
### LOAD DATA
######################################################################################################################################################
### LOAD QUANTIF AND SCORES
QUANTIFandSCORES_LIST = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIFandSCORES_LIST.RDS"))

### LOAD PROFILE MATRIX
PROFMAT_H3K36me2_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/H3K36me2_chipseq_filt_sort_RPGC_profmat.RDS"))
PROFMAT_H3K36me3_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/H3K36me3_chipseq_filt_sort_RPGC_profmat.RDS"))
PROFMAT_K27C_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/K27C_trimmed_filt_sort_RPGC_profmat.RDS"))

######################################################################################################################################################
### PLOT HEATMAP
######################################################################################################################################################

#### FIGURE_1_D_HEATMAP_ranked_Q_H3K36me3_RPGC_GB

rangeheatmap = c(1:1000)
pdf(paste0(workdir, "FIGURES/FIGURE_1_D_E_HEATMAP_PROFILE_CHIPSEQ/FIGURE_1_D_HEATMAP_ranked_Q_H3K36me3_RPGC_GB.pdf"), width=4, height=8)

heatMatrixMat(PROFMAT_H3K36me3_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_H3K36me3_RPGC_GB),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_H3K36me3_RPGC_GB", legend.name = "Q_H3K36me3_RPGC_GB")
heatMatrixMat(PROFMAT_H3K36me2_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_H3K36me3_RPGC_GB),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_H3K36me2_RPGC_GB", legend.name = "Q_H3K36me3_RPGC_GB")
heatMatrixMat(PROFMAT_K27C_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_H3K36me3_RPGC_GB),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_K27C_RPGC_GB", legend.name = "Q_H3K36me3_RPGC_GB")

dev.off()


#### FIGURE_1_E_HEATMAP_ranked_Q_H3K36me2_RPGC_GB

rangeheatmap = c(1:1000)
pdf(paste0(workdir, "FIGURES/FIGURE_1_D_E_HEATMAP_PROFILE_CHIPSEQ/FIGURE_1_E_HEATMAP_ranked_Q_H3K36me2_RPGC_GB.pdf"), width=4, height=8)

heatMatrixMat(PROFMAT_H3K36me3_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_H3K36me2_RPGC_GB),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_H3K36me3_RPGC_GB", legend.name = "Q_H3K36me2_RPGC_GB")
heatMatrixMat(PROFMAT_H3K36me2_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_H3K36me2_RPGC_GB),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_H3K36me2_RPGC_GB", legend.name = "Q_H3K36me2_RPGC_GB")
heatMatrixMat(PROFMAT_K27C_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_H3K36me2_RPGC_GB),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_K27C_RPGC_GB", legend.name = "Q_H3K36me2_RPGC_GB")

dev.off()

# end
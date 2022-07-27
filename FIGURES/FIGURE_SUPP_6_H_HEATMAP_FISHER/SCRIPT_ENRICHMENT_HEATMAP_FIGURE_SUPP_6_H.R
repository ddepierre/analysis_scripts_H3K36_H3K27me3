######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################


# SCRIPT PLOT ENRICHMENT HEATMAP FIGURE 6 B C
# Heatmap of fisher exact test enrichment between genes with specific increase of H3K27me3 in Mes-4 KD or HypB KD (direct Zscore between Mes-4 and HypB condition) and 
# either presence or not of insulator (-> Increase of H3K27me3 detected in Mes-4 KD specifally are associated with the presence of CP190)
# either increase of H3K27me3 detected in interaction mutant Beaf-CP190 (-> Increase of H3K27me3 detected in Mes-4 KD specifally recapitulates increase of H3K27me3 observed in a condition where insulator Beaf is not able to interact with CP190 anymore)


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
### LOAD DATA
######################################################################################################################################################
### LOAD QUANTIF AND SCORES
QUANTIFandSCORES_LIST = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIFandSCORES_LIST.RDS"))
str(QUANTIFandSCORES_LIST)


######################################################################################################################################################
### LOAD FUNCTION
######################################################################################################################################################
# Function to compute fisher exact test from 2 List of vector containing IDs and the plot results of fisher exact test as heatmaps
source(paste0(workdir, "FIGURES/R_PLOT_FUNCTION/FisherList_Hmap.R"))


######################################################################################################################################################
### 2/ FIGURE 6 C
######################################################################################################################################################
### PREPARE DATA TO PLOT
######################################################################################################################################################

ZSCORE_K27M_K27H = QUANTIFandSCORES_LIST$ZSCORE_K27M_K27H
ZSCORE_H3K27me3_muBeaf_exp2 = QUANTIFandSCORES_LIST$ZSCORE_H3K27me3_muBeaf_exp2

Nsplit = 5
ZSCORE_K27M_K27H_GN_splitted = split(names(ZSCORE_K27M_K27H), ceiling(seq_along(names(ZSCORE_K27M_K27H))/ceiling(length(names(ZSCORE_K27M_K27H))/Nsplit)))
ZSCORE_K27M_K27H_VALUE_splitted = split(ZSCORE_K27M_K27H, ceiling(seq_along(ZSCORE_K27M_K27H)/ceiling(length(ZSCORE_K27M_K27H)/Nsplit)))

Nsplit = 5
ZSCORE_H3K27me3_muBeaf_exp2_GN_splitted = split(names(ZSCORE_H3K27me3_muBeaf_exp2), ceiling(seq_along(names(ZSCORE_H3K27me3_muBeaf_exp2))/ceiling(length(names(ZSCORE_H3K27me3_muBeaf_exp2))/Nsplit)))
ZSCORE_H3K27me3_muBeaf_exp2_VALUE_splitted = split(ZSCORE_H3K27me3_muBeaf_exp2, ceiling(seq_along(ZSCORE_H3K27me3_muBeaf_exp2)/ceiling(length(ZSCORE_H3K27me3_muBeaf_exp2)/Nsplit)))

# compute fisher exact test

ZSCORE_K27M_K27H_ZSCORE_H3K27me3_muBeaf_exp2 = fisher_namesList(ZSCORE_K27M_K27H_GN_splitted, ZSCORE_H3K27me3_muBeaf_exp2_GN_splitted, length(ZSCORE_K27M_K27H))


######################################################################################################################################################
### PLOT MATRICES
######################################################################################################################################################
  TH=1e-10 # pvalue vizualisation threshold

Hmap_pval(ZSCORE_K27M_K27H_ZSCORE_H3K27me3_muBeaf_exp2, values_ordonnée = ZSCORE_K27M_K27H_VALUE_splitted, values_abscisse = ZSCORE_H3K27me3_muBeaf_exp2_VALUE_splitted,
						paste0(workdir, "FIGURES/FIGURE_6_B_C_HEATMAP_FISHER/HEATMAP_ENRICHMENT_FISHER_ZSCORE_K27M_K27H_ZSCORE_H3K27me3_muBeaf_exp2_FIGURE_SUPP_6_H_", length(ZSCORE_K27M_K27H),".pdf"),T,threshold=TH, odds_scale = c(0,2))







#end

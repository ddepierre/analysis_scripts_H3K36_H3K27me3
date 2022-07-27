######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################


# SCRIPT PLOT ENRICHMENT HEATMAP FIGURE 3 C
# Heatmap of fisher exact test enrichment between de-regulated genes and genes with increase of H3K27me3
# Statistical validation that Mes-4 and HypB protect genes from down-regulation by a H3K27me3 spreading.

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
library("SuperExactTest")
library("Vennerable")

## define workdir
workdir = "/path/to/analysis_scripts_H3K36_H3K27me3/" 


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
### PREPARE DATA TO PLOT
######################################################################################################################################################

commonGeneNames_HypB = Reduce(intersect,list(names(QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C), names(QUANTIFandSCORES_LIST$LOGFC_HYPBKD_RPGC_pval005)))
commonGeneNames_Mes4 = Reduce(intersect,list(names(QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C), names(QUANTIFandSCORES_LIST$LOGFC_MES4KD_RPGC_pval005)))

ZSCORE_K27H_K27C_f = QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C[names(QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C) %in% commonGeneNames_HypB]
ZSCORE_K27M_K27C_f = QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C[names(QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C) %in% commonGeneNames_Mes4]

DIFFEXPR_HYPB = QUANTIFandSCORES_LIST$LOGFC_HYPBKD_RPGC_pval005[names(QUANTIFandSCORES_LIST$LOGFC_HYPBKD_RPGC_pval005) %in% commonGeneNames_HypB]
DIFFEXPR_MES4 = QUANTIFandSCORES_LIST$LOGFC_MES4KD_RPGC_pval005[names(QUANTIFandSCORES_LIST$LOGFC_MES4KD_RPGC_pval005) %in% commonGeneNames_Mes4]

################################################################################################################################################################################################
### PLOT MATRICES
################################################################################################################################################################################################
ZSCORE_K27H_K27C_f
ZSCORE_K27M_K27C_f

Nsplit = 5
ZSCORE_K27H_K27C_GN_splitted = split(names(ZSCORE_K27H_K27C_f), ceiling(seq_along(names(ZSCORE_K27H_K27C_f))/ceiling(length(names(ZSCORE_K27H_K27C_f))/Nsplit)))
ZSCORE_K27H_K27C_VALUE_splitted = split(ZSCORE_K27H_K27C_f, ceiling(seq_along(ZSCORE_K27H_K27C_f)/ceiling(length(ZSCORE_K27H_K27C_f)/Nsplit)))

Nsplit = 5
ZSCORE_K27M_K27C_GN_splitted = split(names(ZSCORE_K27M_K27C_f), ceiling(seq_along(names(ZSCORE_K27M_K27C_f))/ceiling(length(names(ZSCORE_K27M_K27C_f))/Nsplit)))
ZSCORE_K27M_K27C_VALUE_splitted = split(ZSCORE_K27M_K27C_f, ceiling(seq_along(ZSCORE_K27M_K27C_f)/ceiling(length(ZSCORE_K27M_K27C_f)/Nsplit)))

Nsplit = 5
DIFFEXPR_HYPB_GN_splitted = split(names(DIFFEXPR_HYPB), ceiling(seq_along(names(DIFFEXPR_HYPB))/ceiling(length(names(DIFFEXPR_HYPB))/Nsplit)))
DIFFEXPR_HYPB_VALUE_splitted = split(DIFFEXPR_HYPB, ceiling(seq_along(DIFFEXPR_HYPB)/ceiling(length(DIFFEXPR_HYPB)/Nsplit)))

Nsplit = 5
DIFFEXPR_MES4_GN_splitted = split(names(DIFFEXPR_MES4), ceiling(seq_along(names(DIFFEXPR_MES4))/ceiling(length(names(DIFFEXPR_MES4))/Nsplit)))
DIFFEXPR_MES4_VALUE_splitted = split(DIFFEXPR_MES4, ceiling(seq_along(DIFFEXPR_MES4)/ceiling(length(DIFFEXPR_MES4)/Nsplit)))



ZSCORE_K27H_K27C__DIFFEXPR_HYPB = fisher_namesList(ZSCORE_K27H_K27C_GN_splitted, DIFFEXPR_HYPB_GN_splitted, length(ZSCORE_K27H_K27C_f))
TH=1e-10
Hmap_pval(ZSCORE_K27H_K27C__DIFFEXPR_HYPB, values_ordonnée = ZSCORE_K27H_K27C_VALUE_splitted, values_abscisse = DIFFEXPR_HYPB_VALUE_splitted,
						paste0(workdir, "FIGURES/FIGURE_3_C_HEATMAP_FISHER/HEATMAP_ENRICHMENT_FISHER_ZSCORE_K27H_K27C__DIFFEXPR_HYPB_FIGURE_3_C", length(ZSCORE_K27H_K27C_f),".pdf"),T,threshold=TH, odds_scale = c(0,2))


ZSCORE_K27M_K27C__DIFFEXPR_MES4 = fisher_namesList(ZSCORE_K27M_K27C_GN_splitted, DIFFEXPR_MES4_GN_splitted, length(ZSCORE_K27M_K27C_f))
TH=1e-10
Hmap_pval(ZSCORE_K27M_K27C__DIFFEXPR_MES4, values_ordonnée = ZSCORE_K27M_K27C_VALUE_splitted, values_abscisse = DIFFEXPR_MES4_VALUE_splitted,
						paste0(workdir, "FIGURES/FIGURE_3_C_HEATMAP_FISHER/HEATMAP_ENRICHMENT_FISHER_ZSCORE_K27M_K27C__DIFFEXPR_MES4_FIGURE_3_C", length(ZSCORE_K27M_K27C_f),".pdf"),T,threshold=TH, odds_scale = c(0,2))






######################################################################################################################################################
### END
######################################################################################################################################################

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
# Heatmap of fisher exact test enrichment between genes with increase of H3K27me3 and genes classified accroding to their position to H3K27me3 domains border.
# Statistical validation that genes silenced by a H3K27me3 spreading when H3K36 HMT are depleted are localised in a close neighborhood of H3K27me3 domains border.

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
workdir = ""


######################################################################################################################################################
### LOAD DATA
######################################################################################################################################################
### LOAD QUANTIF AND SCORES
QUANTIFandSCORES_LIST = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIFandSCORES_LIST.RDS"))
str(QUANTIFandSCORES_LIST)

### LOAD GENES ID LIST CLASSIFIED ON THEIR POSITION TO H3K27ME3 DOMAINS BORDER (done in DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/NORMR_H3K27me3_DOMAINS.R)
List_genes_DOM = readRDS(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/List_genes_DOM.RDS"))
str(List_genes_DOM)

# Load ref genes Granges file
gene_dm6_gr = readRDS(paste0(workdir, "DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.RDS"))
names(gene_dm6_gr) = paste0(names(gene_dm6_gr), ".1")

######################################################################################################################################################
### LOAD FUNCTION
######################################################################################################################################################
# Function to compute fisher exact test from 2 List of vector containing IDs and the plot results of fisher exact test as heatmaps
source(paste0(workdir, "FIGURES/R_PLOT_FUNCTION/FisherList_Hmap.R"))

## extact a vector of name corresponding to X% of an ordered value
getNameList = function(Vec, topdown = "top", prct = 10){
  Vec = Vec[order(Vec, decreasing=T)]
  if(topdown %in% "top"){
    GN = names(Vec[Vec > quantile(Vec, (100-prct)/100)])
  }
  if(topdown %in% "down"){
     GN = names(Vec[Vec < quantile(Vec, (prct)/100)])
  }
  if(topdown %in% "mid"){
     tmp1 = names(Vec[Vec < quantile(Vec, (100/2-prct/2)/100)])
     tmp2 = names(Vec[Vec < quantile(Vec, (100/2-prct/2+prct)/100)])
     GN = tmp2[tmp2 %ni% tmp1]
  }
  return(GN)
}


######################################################################################################################################################
### PREPARE DATA TO PLOT
######################################################################################################################################################


# FILTER GENES WITH ACTIVE GENS
str(List_genes_DOM)
List_genes_DOM_f = lapply(List_genes_DOM, function(feat1){feat1 = feat1[feat1 %in% QUANTIFandSCORES_LIST$ACTIVE_GENES]})
str(List_genes_DOM_f)

# get ZSCORE H3K27me3 vector for active genes
ZSCORE_K27H_K27C_f = QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C[names(QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES]
ZSCORE_K27M_K27C_f = QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C[names(QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES]

# get H3K27me3 increase genes as Zscore top 20%
ZSCORE_K27H_K27C_f_top20 = getNameList(ZSCORE_K27H_K27C_f, topdown = "top", prct = 20 )
ZSCORE_K27M_K27C_f_top20 = getNameList(ZSCORE_K27M_K27C_f, topdown = "top", prct = 20 )

# Create a list of Genes ID with 4 classes as:
LIST_UP_K27_SHARED_GN = list(
  HM = ZSCORE_K27H_K27C_f_top20[ZSCORE_K27H_K27C_f_top20 %in% ZSCORE_K27M_K27C_f_top20], # increase of H3K27me3 in both Mes-4 KD and HypB KD
  H = ZSCORE_K27H_K27C_f_top20[ZSCORE_K27H_K27C_f_top20 %ni% ZSCORE_K27M_K27C_f_top20], # increase of H3K27me3 only in HypB KD
  M = ZSCORE_K27M_K27C_f_top20[ZSCORE_K27M_K27C_f_top20 %ni% ZSCORE_K27H_K27C_f_top20], # increase of H3K27me3 only in Mes-4 KD
  random = sample(names(ZSCORE_K27M_K27C_f), 700) # Random set of genes as control
)

# compute fisher exact test
List_genesDOM_LIST_UP_K27_SHARED_GN = fisher_namesList(list1 = List_genes_DOM_f, list2 = LIST_UP_K27_SHARED_GN, total = length(unlist(List_genes_DOM_f)), test.side = "greater")
# result of fisher exact test can be visualized as raw matrix : 


######################################################################################################################################################
### PLOT MATRICES
######################################################################################################################################################
  TH=1e-10 # pvalue vizualisation threshold
  Hmap_pval(List_genesDOM_LIST_UP_K27_SHARED_GN, values_ordonnée = List_genes_DOM_f, values_abscisse = LIST_UP_K27_SHARED_GN, odds_scale = c(0,3),
						paste0(workdir, "/FIGURES/FIGURE_3_C_HEATMAP_FISHER/HEATMAP_ENRICHMENT_FISHER_UPK27_vs_GENESpositionToBorders_FIGURE_4_C.pdf"),T,threshold=TH)




























#end

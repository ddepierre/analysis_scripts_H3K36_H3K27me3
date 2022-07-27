######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# SCRIPT PLOT VENN DIAGRAM FIGURE 3 A
# Description of 3 classes of genes that show a H3K27me3 increase upon H3K36 HMT depletion
# Genes specifically under HypB (H) or Mes-4 (M) regulation or genes regulated the same by both depletion (HM)

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

# Load ref genes Granges file
gene_dm6_gr = readRDS(paste0(workdir, "DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.RDS"))
names(gene_dm6_gr) = paste0(names(gene_dm6_gr), ".1")

######################################################################################################################################################
### LOAD FUNCTION
######################################################################################################################################################
# Function to plot venn diagram using Vennerable and SuperExactTest libraries
source(paste0(workdir, "FIGURES/R_PLOT_FUNCTION/VENN_vennerable.R"))

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
### PREPARE GENES ID LIST TO PLOT
######################################################################################################################################################

## Get matrice of overlapping genes (by coordinates)
OVLP_GENES = findOverlaps(gene_dm6_gr,gene_dm6_gr)
OVLP_GENES_MAT = cbind(OVLP_GENES@from, OVLP_GENES@to) # get matrix with genes ID
dim(OVLP_GENES_MAT)
OVLP_GENES_MAT = OVLP_GENES_MAT[OVLP_GENES_MAT[,1] != OVLP_GENES_MAT[,2],] # keep only 1 time the gene id if it overlaps
dim(OVLP_GENES_MAT)


## From genes that overlap with others, keep only 1 gene over the X overlapping
OVLP_GENES_MAT_tokeep = OVLP_GENES_MAT
GENES_OVLP_tokeep = c()
for(i in 1:length(gene_dm6_gr)){
  # print(i)
  # print(OVLP_GENES_MAT_tokeep[which(OVLP_GENES_MAT_tokeep[,1] %in% i + OVLP_GENES_MAT_tokeep[,2] %in% i == 1),])

  if(is.na(unique(c(OVLP_GENES_MAT_tokeep[which(OVLP_GENES_MAT_tokeep[,1] %in% i + OVLP_GENES_MAT_tokeep[,2] %in% i == 1),]))[1]) %in% F){
    GENES_OVLP_tokeep = c(GENES_OVLP_tokeep, unique(c(OVLP_GENES_MAT_tokeep[which(OVLP_GENES_MAT_tokeep[,1] %in% i + OVLP_GENES_MAT_tokeep[,2] %in% i == 1),]))[1])
    OVLP_GENES_MAT_tokeep = OVLP_GENES_MAT_tokeep[-which(OVLP_GENES_MAT_tokeep[,1] %in% i + OVLP_GENES_MAT_tokeep[,2] %in% i == 1),]
  }
}
# GENES_OVLP_tokeep is the vector of genes ID that overlaps but only 1 per overlapping hub is kept

gene_dm6_gr_NOOVLP = gene_dm6_gr[-unique(OVLP_GENES_MAT[,1])] # remove genes that overlap = get genes with no overlap
gene_dm6_gr_NOOVLPtoadd = gene_dm6_gr[GENES_OVLP_tokeep] # get genes that overlap where we kept only 1 gene per "overlapping hub"
gene_dm6_gr_NOOVLP = c(gene_dm6_gr_NOOVLP, gene_dm6_gr_NOOVLPtoadd) # combine genes with no overlap with other + 1 gene per overlapping hub

# gene_dm6_gr_NOOVLP are unique gene sites, i.e. if 3 genes are overlapping, only 1 position is kept over the 3 not to biase the enrichment test


ZSCORE_K27H_K27C_f = QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C[names(QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES] # get ZSCORE H3K27me3 vector for active genes
ZSCORE_K27M_K27C_f = QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C[names(QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES] # get ZSCORE H3K27me3 vector for active genes

UP_ZSCORE_K27H_K27C_f = getNameList(ZSCORE_K27H_K27C_f, "top", 20) # get H3K27me3 increase genes as top 20%
UP_ZSCORE_K27M_K27C_f = getNameList(ZSCORE_K27M_K27C_f, "top", 20) # get H3K27me3 increase genes as top 20%

UP_ZSCORE_K27H_K27C_f_noOvlp = UP_ZSCORE_K27H_K27C_f[UP_ZSCORE_K27H_K27C_f %in% names(gene_dm6_gr_NOOVLP)]# keep only no ovlp genes
UP_ZSCORE_K27M_K27C_f_noOvlp = UP_ZSCORE_K27M_K27C_f[UP_ZSCORE_K27M_K27C_f %in% names(gene_dm6_gr_NOOVLP)]# keep only no ovlp genes

# put the 2 vectors as a list
List1 = list(
              UP_ZSCORE_K27H_K27C_f_noOvlp = UP_ZSCORE_K27H_K27C_f_noOvlp,
              UP_ZSCORE_K27M_K27C_f_noOvlp = UP_ZSCORE_K27M_K27C_f_noOvlp
            )


######################################################################################################################################################
### PLOT VENN
######################################################################################################################################################
outfig = paste0(workdir, "FIGURES/FIGURE_4_B_VENN_DIAGRAM/")

plotVENNerable(data_to_venn = List1, Nb_REF = length(QUANTIFandSCORES_LIST$ACTIVE_GENES), outdir = outfig)


#end

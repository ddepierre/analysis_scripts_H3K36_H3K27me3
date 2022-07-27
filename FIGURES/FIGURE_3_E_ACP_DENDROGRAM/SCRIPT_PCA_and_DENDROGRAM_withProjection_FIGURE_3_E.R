######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# SCRIPT TO DO AND PLOT PCA AND DENDROGRAM FROM PCA PROJETION
#  This analysis and visualization sum up the global relation between scores and quantifications

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
workdir = "/home/depierre/Bureau/analysis_scripts_H3K36_H3K27me3/"


######################################################################################################################################################
### LOAD DATA
######################################################################################################################################################
### LOAD QUANTIF AND SCORES
QUANTIFandSCORES_LIST = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIFandSCORES_LIST.RDS"))
str(QUANTIFandSCORES_LIST)


# Load ref genes Granges file
gene_dm6_gr = readRDS(paste0(workdir, "DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.RDS"))
names(gene_dm6_gr) = paste0(names(gene_dm6_gr), ".1")
gene_dm6_gr

######################################################################################################################################################
### LOAD FUNCTION
######################################################################################################################################################
# Function to compute PCA for list of numeric vectors
source(paste0(workdir, "FIGURES/R_PLOT_FUNCTION/PCA_and_DENDROGRAM.R"))

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


###################################################################################################################################################
###  PREPARE LIST OF DATA FOR PCA
###################################################################################################################################################

LIST_DOWNREG = list(
	LOG2RATIO_K36me3_me2_f = QUANTIFandSCORES_LIST$LOG2RATIO_K36me3_me2[names(QUANTIFandSCORES_LIST$LOG2RATIO_K36me3_me2) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
	LOG2RATIO_K36me2_me3_f = QUANTIFandSCORES_LIST$LOG2RATIO_K36me3_me2[names(QUANTIFandSCORES_LIST$LOG2RATIO_K36me3_me2) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES]*-1,
	Q_H3K36me3_RPGC_GB_f = QUANTIFandSCORES_LIST$Q_H3K36me3_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_H3K36me3_RPGC_GB) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
	Q_H3K36me2_RPGC_GB_f = QUANTIFandSCORES_LIST$Q_H3K36me2_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_H3K36me2_RPGC_GB) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
	ZSCORE_K27H_K27C_f = QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C[names(QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
	ZSCORE_K27M_K27C_f = QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C[names(QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
  Q_H3K27me3_GB_f = QUANTIFandSCORES_LIST$Q_K27C_RPGC_GB[names(QUANTIFandSCORES_LIST$Q_K27C_RPGC_GB) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
	DNREG_HYPBKD = QUANTIFandSCORES_LIST$LOGFC_HYPBKD_RPGC_pval005[names(QUANTIFandSCORES_LIST$LOGFC_HYPBKD_RPGC_pval005) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES]*-1, # invert logFC sign to have the highest logFC as the Down-regulated
	DNREG_MES4KD = QUANTIFandSCORES_LIST$LOGFC_MES4KD_RPGC_pval005[names(QUANTIFandSCORES_LIST$LOGFC_MES4KD_RPGC_pval005) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES]*-1, # invert logFC sign to have the highest logFC as the Down-regulated
	dist_gene2border_f = QUANTIFandSCORES_LIST$dist_gene2border[names(QUANTIFandSCORES_LIST$dist_gene2border) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES]*-1,
  Q_H3K27ac_rep1_GB_f = QUANTIFandSCORES_LIST$Q_H3K27ac_rep1_GB[names(QUANTIFandSCORES_LIST$Q_H3K27ac_rep1_GB) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
  Q_H3K4me1_rep1_GB_f = QUANTIFandSCORES_LIST$Q_H3K4me1_rep1_GB[names(QUANTIFandSCORES_LIST$Q_H3K4me1_rep1_GB) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
  Q_H3K4me3_rep1_GB_f = QUANTIFandSCORES_LIST$Q_H3K4me3_rep1_GB[names(QUANTIFandSCORES_LIST$Q_H3K4me3_rep1_GB) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
  Q_H3_Rep1 = QUANTIFandSCORES_LIST$Q_eGFPdsRNA_H3_ip_Rep1_GB[names(QUANTIFandSCORES_LIST$Q_eGFPdsRNA_H3_ip_Rep1_GB) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES],
  Q_H3K9Me3 = QUANTIFandSCORES_LIST$Q_GFPRNAi_H3K9Me3ChIP3_GB[names(QUANTIFandSCORES_LIST$Q_GFPRNAi_H3K9Me3ChIP3_GB) %in% QUANTIFandSCORES_LIST$ACTIVE_GENES]
)



###################################################################################################################################################
###  DO AND PLOT PCA
###################################################################################################################################################


pdf(paste0(workdir, "FIGURES/FIGURE_3_E_ACP_DENDROGRAM/ACP_DOWNREG_genes_FIGURE_3_E.pdf"), height = 10, width= 10)
do_acp_and_clustering(LIST_DOWNREG,
                      list.clust=list(c(1:3)),
                      list.cpa=list(c(1,2),c(1,6)),
                      cpa.title="PCA", dist.method="euclidean",
                      choix="var",
                      col.ind = "1")
dev.off()




#end

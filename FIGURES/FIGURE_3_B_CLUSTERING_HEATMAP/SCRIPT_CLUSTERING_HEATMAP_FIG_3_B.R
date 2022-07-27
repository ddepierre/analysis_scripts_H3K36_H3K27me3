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
### prepare matrix of Diff expr analysis logFC for heatmap and clustering
######################################################################################################################################################


LOGFC_DiffExpr_HYPBKD_RPGC_pval005 = as.matrix(QUANTIFandSCORES_LIST$LOGFC_HYPBKD_RPGC_pval005)
LOGFC_DiffExpr_MES4KD_RPGC_pval005 = as.matrix(QUANTIFandSCORES_LIST$LOGFC_MES4KD_RPGC_pval005)
comRownames = Reduce(intersect,list(rownames(LOGFC_DiffExpr_HYPBKD_RPGC_pval005), rownames(LOGFC_DiffExpr_MES4KD_RPGC_pval005)))
LOGFC_DiffExpr_HYPBKD_RPGC_pval005 = LOGFC_DiffExpr_HYPBKD_RPGC_pval005[rownames(LOGFC_DiffExpr_HYPBKD_RPGC_pval005) %in% comRownames,1,drop=F]
LOGFC_DiffExpr_MES4KD_RPGC_pval005 = LOGFC_DiffExpr_MES4KD_RPGC_pval005[rownames(LOGFC_DiffExpr_MES4KD_RPGC_pval005) %in% comRownames,1,drop=F]
LOGFC_DiffExpr_HYPBKD_RPGC_pval005 = LOGFC_DiffExpr_HYPBKD_RPGC_pval005[comRownames, 1,drop=F]
LOGFC_DiffExpr_MES4KD_RPGC_pval005 = LOGFC_DiffExpr_MES4KD_RPGC_pval005[comRownames, 1,drop=F]
mat_LOGFC_RPGC_DE_pval005 = matrix(c(LOGFC_DiffExpr_MES4KD_RPGC_pval005[,1], LOGFC_DiffExpr_HYPBKD_RPGC_pval005[,1]), ncol=2)
rownames(mat_LOGFC_RPGC_DE_pval005) = comRownames
colnames(mat_LOGFC_RPGC_DE_pval005) = c("DE_MES4", "DE_HYPB")


###################################################################################################################################################
###  Plot heatmap with clustering of genes according to their differential expression in Mes-4 and HypB
###################################################################################################################################################

## 1/ FILTER mat_ZSCORE_DE matrix
mat = mat_LOGFC_RPGC_DE_pval005
TH = 0.3 # threshold du LOGFC
GNnames_filt = unique(c(rownames(mat[abs(c(mat[,1])) > TH,]), rownames(mat[abs(c(mat[,2])) > TH,])))
mat_f = mat[rownames(mat) %in% GNnames_filt,,drop=F]
mat_LOGFC_RPGC_DE_pval005_f = mat_f

## 2/ CUT OUTLIERS
mat = mat_LOGFC_RPGC_DE_pval005_f
TH = 0.02 ## let 2% outliers be erase
mat_temp = mat
mat[mat[,1] < quantile(c(mat_temp), TH),1] = quantile(c(mat_temp), TH)
mat[mat[,2] < quantile(c(mat_temp), TH),2] = quantile(c(mat_temp), TH)
mat[mat[,1] > quantile(c(mat_temp), 1-TH),1] = quantile(c(mat_temp), 1-TH)
mat[mat[,2] > quantile(c(mat_temp), 1-TH),2] = quantile(c(mat_temp), 1-TH)
mat_LOGFC_RPGC_DE_pval005_f_rmol = mat

## 3/ CLUSTER GENES AND PLOT HEATMAP
## Row- and column-wise clustering
hr <- hclust(dist(mat_LOGFC_RPGC_DE_pval005_f_rmol, method = "euclidean", diag = FALSE, upper = FALSE), method="complete")
hc <- hclust(as.dist(1-cor(mat_LOGFC_RPGC_DE_pval005_f_rmol, method="pearson")), method="complete")
## Tree cutting
# mycl <- cutree(hr, h=max(hr$height)/1.5)
mycl <- cutree(hr, h=max(hr$height)/2.3)
mycolhr <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
mycolhr <- mycolhr[as.vector(mycl)]
## Plot heatmap
mycol <- colorpanel(100, "red4", "black", "green3") # or try redgreen(75)
pdf(paste0(workdir, "FIGURES/FIGURE_3_B_CLUSTERING_HEATMAP/HEATMAP_CLUSTERING_DEA_mat_LOGFC_RPGC_DE_pval005_FIGURE_3_B.pdf"))
# heatmap.2(mat_LOGFC_RPGC_DE_pval005_f_rmol, col=mycol,density.info="density", trace="none", RowSideColors=mycolhr, breaks = seq(quantile(c(mat_temp), TH), quantile(c(mat_temp), 1-TH), length.out = 101))
# heatmap.2(mat_LOGFC_RPGC_DE_pval005_f_rmol, col=mycol,density.info="density", trace="none", RowSideColors=mycolhr, breaks = seq(-10,10, length.out = 101))
heatmap.2(mat_LOGFC_RPGC_DE_pval005_f_rmol, col=mycol,density.info="non", trace="none", RowSideColors=mycolhr)
heatmap.2(mat_LOGFC_RPGC_DE_pval005_f_rmol, col=mycol, scale = "col",density.info="density", trace="none", RowSideColors=mycolhr)
dev.off()


###################################################################################################################################################
###  Export and annotate clustered genes
###################################################################################################################################################
# get a list of no up nor down-reg for control 
noUP_noDN = QUANTIFandSCORES_LIST$ACTIVE_GENES[QUANTIFandSCORES_LIST$ACTIVE_GENES %ni% rownames(mat_LOGFC_RPGC_DE_pval005)]


DEA_Cluster_7G = list(
noUP_noDN = sample(noUP_noDN, 1036),
U7 = names(mycl[mycl %in% 1]),
U6 = names(mycl[mycl %in% 4]),
U5 = names(mycl[mycl %in% 3]),
D4 = names(mycl[mycl %in% 2]),
D3 = names(mycl[mycl %in% 6]),
D2 = names(mycl[mycl %in% 7]),
D1 = names(mycl[mycl %in% 5]))

saveRDS(DEA_Cluster_7G, paste0(workdir, "FIGURES/FIGURE_3_B_CLUSTERING_HEATMAP/CLUSTERING_DEA_genes.RDS"))

###################################################################################################################################################
###  END
###################################################################################################################################################

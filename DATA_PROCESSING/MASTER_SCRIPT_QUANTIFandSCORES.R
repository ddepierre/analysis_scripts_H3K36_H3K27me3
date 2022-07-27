######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# MASTER SCRIPT QUANTIF & SCORES
## This Script generates a list saved in RDS format (R object data) an containing all the quantifications and scores used in the figures.
## See DATA_PROCESSING/ folder for scripts used to obtain these quantif
## Quantifications from other published studies are obtained the exact same way but scripts are not detailled in the DATA_PROCESSING/ folder

workdir = "/home/depierre/Bureau/analysis_scripts_H3K36_H3K27me3/" 

QUANTIFandSCORES_LIST = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIFandSCORES_LIST.RDS"))




######################################################################################################################################################
######################################################################################################################################################
# R version 3.4.2 (2017-09-28) -- "Short Summer"
# Copyright (C) 2017 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

######################################################################################################################################################
### LIBRARY
######################################################################################################################################################
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")


#########################################################################################################################
### FUNCTION
#########################################################################################################################
'%ni%' = Negate('%in%')

#########  FUNCTION
### -> Compute Zscore as difference normalized by sqrt of means
matMeans <- function(X,Y){ mean(c(X,Y)) }

matReplaceNA = function(M){
  if(is.null(dim(M))){
    Mv = M
    Mv[which(is.na(Mv))] = 0
    Mv[is.finite(Mv) %in% F] = 0
  }else{
    Mv = c(M)
    Mv[which(is.na(Mv))] = 0
    Mv[is.finite(Mv) %in% F] = 0
    Mv = matrix(Mv, ncol=1)
    rownames(Mv) = rownames(M)
  }
  return(Mv)
}


computeZscore = function(Q_KD, Q_CTRL){
    comRownames = Reduce(intersect,list(names(Q_KD), names(Q_CTRL)))
    Q_KD = Q_KD[names(Q_KD) %in% comRownames]
    Q_CTRL = Q_CTRL[names(Q_CTRL) %in% comRownames]
    ZSCORE = (Q_KD-Q_CTRL[names(Q_KD)])/sqrt(matrix(mapply(matMeans, Q_KD, Q_CTRL), ncol=1))
    ZSCORE = matReplaceNA(ZSCORE)
    rownames(ZSCORE) = names(Q_KD)
    ZSCORE  = ZSCORE[,1]
    ZSCORE = ZSCORE[order(ZSCORE, decreasing=T),drop=F]
    return(ZSCORE)
}


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


#########################################################################################################################
# INITIATE LIST
#########################################################################################################################
QUANTIFandSCORES_LIST = list()

#########################################################################################################################
# get genes ref from UCSC
#########################################################################################################################
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
gene_dm6_gr <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
names(gene_dm6_gr) = paste0(names(gene_dm6_gr),".1")
gene_dm6_gr <- gene_dm6_gr[seqnames(gene_dm6_gr) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4",  "chrX",  "chrY")]
gene_dm6_gr_TSS = resize(gene_dm6_gr, 1, "start")

#########################################################################################################################
## QUANTIF GENE BODY DONE IN DATA_PROCESSING/QUANTIF/GENE_BODY/
#########################################################################################################################
workdir = "/home/depierre/Bureau/analysis_scripts_H3K36_H3K27me3/" # on CUVIER10

##########################################################
####   CHIP SEQ K27     
# LOAD QUANTIF
Q_C27_1_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/C27_1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H27_2_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/H27_2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_K27C_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/K27C_RPGC_readsCounts_GB_SCALED.RDS"))
Q_K27H_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/K27H_RPGC_readsCounts_GB_SCALED.RDS"))
Q_K27M_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/K27M_RPGC_readsCounts_GB_SCALED.RDS"))
Q_M27_2_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/M27_2_RPGC_readsCounts_GB_SCALED.RDS"))


## COMPUTE ZSCORE
ZSCORE_H27_2_C27_1 = computeZscore(Q_H27_2_RPGC_GB, Q_C27_1_RPGC_GB)
ZSCORE_K27H_K27C = computeZscore(Q_K27H_RPGC_GB, Q_K27C_RPGC_GB)
ZSCORE_M27_2_C27_1 = computeZscore(Q_M27_2_RPGC_GB, Q_C27_1_RPGC_GB)
ZSCORE_K27M_K27C = computeZscore(Q_K27M_RPGC_GB, Q_K27C_RPGC_GB)
# direct differential score between Mes-4 KD and HypB KD to extract specific H3K27me3 increases
ZSCORE_K27M_K27H = computeZscore(Q_K27M_RPGC_GB, Q_K27H_RPGC_GB)


# ORDER QUANTIF decreasingly
Q_C27_1_RPGC_GB = Q_C27_1_RPGC_GB[order(Q_C27_1_RPGC_GB, decreasing=T),drop=F]
Q_K27C_RPGC_GB = Q_K27C_RPGC_GB[order(Q_K27C_RPGC_GB, decreasing=T),drop=F]
Q_H27_2_RPGC_GB = Q_H27_2_RPGC_GB[order(Q_H27_2_RPGC_GB, decreasing=T),drop=F]
Q_K27H_RPGC_GB = Q_K27H_RPGC_GB[order(Q_K27H_RPGC_GB, decreasing=T),drop=F]
Q_M27_2_RPGC_GB = Q_M27_2_RPGC_GB[order(Q_M27_2_RPGC_GB, decreasing=T),drop=F]
Q_K27M_RPGC_GB = Q_K27M_RPGC_GB[order(Q_K27M_RPGC_GB, decreasing=T),drop=F]

# ORDER ZSCORE decreasingly
ZSCORE_K27H_K27C = ZSCORE_K27H_K27C[order(ZSCORE_K27H_K27C, decreasing=T),drop=F]
ZSCORE_H27_2_C27_1 = ZSCORE_H27_2_C27_1[order(ZSCORE_H27_2_C27_1, decreasing=T),drop=F]
ZSCORE_K27M_K27C = ZSCORE_K27M_K27C[order(ZSCORE_K27M_K27C, decreasing=T),drop=F]
ZSCORE_M27_2_C27_1 = ZSCORE_M27_2_C27_1[order(ZSCORE_M27_2_C27_1, decreasing=T),drop=F]

ZSCORE_K27M_K27H = ZSCORE_K27M_K27H[order(ZSCORE_K27M_K27H, decreasing=T),drop=F]


# ADD in QUANTIFandSCORES_LIST
QUANTIFandSCORES_LIST

QUANTIFandSCORES_LIST$Q_C27_1_RPGC_GB = Q_C27_1_RPGC_GB
QUANTIFandSCORES_LIST$Q_K27C_RPGC_GB = Q_K27C_RPGC_GB
QUANTIFandSCORES_LIST$Q_H27_2_RPGC_GB = Q_H27_2_RPGC_GB
QUANTIFandSCORES_LIST$Q_K27H_RPGC_GB = Q_K27H_RPGC_GB
QUANTIFandSCORES_LIST$Q_M27_2_RPGC_GB = Q_M27_2_RPGC_GB
QUANTIFandSCORES_LIST$Q_K27M_RPGC_GB = Q_K27M_RPGC_GB
QUANTIFandSCORES_LIST$ZSCORE_K27H_K27C = ZSCORE_K27H_K27C
QUANTIFandSCORES_LIST$ZSCORE_H27_2_C27_1 = ZSCORE_H27_2_C27_1
QUANTIFandSCORES_LIST$ZSCORE_K27M_K27C = ZSCORE_K27M_K27C
QUANTIFandSCORES_LIST$ZSCORE_M27_2_C27_1 = ZSCORE_M27_2_C27_1
QUANTIFandSCORES_LIST$ZSCORE_K27M_K27H = ZSCORE_K27M_K27H


##########################################################
####   CHIP SEQ K36     
# CHIP SEQ K36 QUANTIF ENTIRE GB SCALED
Q_H3K36me2_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/H3K36me2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me3_RPGC_GB = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/H3K36me3_RPGC_readsCounts_GB_SCALED.RDS"))

# RATIO K36me2/3 3/2
LOG2RATIO_K36me3_me2 = log2(Q_H3K36me3_RPGC_GB+1)-log2(Q_H3K36me2_RPGC_GB[names(Q_H3K36me3_RPGC_GB)]+1)
LOG2RATIO_K36me2_me3 = log2(Q_H3K36me2_RPGC_GB+1)-log2(Q_H3K36me3_RPGC_GB[names(Q_H3K36me2_RPGC_GB)]+1)

# ORDER QUANTIF and ratio decreasingly
Q_H3K36me2_RPGC_GB = Q_H3K36me2_RPGC_GB[order(Q_H3K36me2_RPGC_GB, decreasing=T),drop=F]
Q_H3K36me3_RPGC_GB = Q_H3K36me3_RPGC_GB[order(Q_H3K36me3_RPGC_GB, decreasing=T),drop=F]
LOG2RATIO_K36me3_me2 = LOG2RATIO_K36me3_me2[order(LOG2RATIO_K36me3_me2, decreasing=T),drop=F]
LOG2RATIO_K36me2_me3 = LOG2RATIO_K36me2_me3[order(LOG2RATIO_K36me2_me3, decreasing=T),drop=F]

# ADD in QUANTIFandSCORES_LIST
QUANTIFandSCORES_LIST$Q_H3K36me2_RPGC_GB = Q_H3K36me2_RPGC_GB
QUANTIFandSCORES_LIST$Q_H3K36me3_RPGC_GB = Q_H3K36me3_RPGC_GB
QUANTIFandSCORES_LIST$LOG2RATIO_K36me3_me2 = LOG2RATIO_K36me3_me2
QUANTIFandSCORES_LIST$LOG2RATIO_K36me2_me3 = LOG2RATIO_K36me2_me3


#########################################################################################################################
## PUBLIC DATA QUANTIFICATION (used in FIGURE_3_E_ACP_DENDROGRAM or Supp validation figures)
## FROM FASTQ to GENE BODY QUANTIFICATION, THE SAME PROCESSING HAS BEEN USED for public data and for ChIPseq generated for this project 
#########################################################################################################################

##################################################
## H3K9me3 (GSE99027) and H3 CHIPSEQ (GSE113470)
##################################################
# LOAD QUANTIF
# Q_eGFPdsRNA_H3_ip_Rep1_GB = readRDS("/home/depierre/Bureau/work_cperrois/ALIGN_GSE99027/QUANITF/GB/eGFPdsRNA_H3_ip_Rep1_readsCounts_GB_SCALED.RDS")
# Q_GFPRNAi_H3K9Me3ChIP3_GB = readRDS("/home/depierre/Bureau/work_cperrois/ALIGN_GSE99027/QUANITF/GB/GFPRNAi_H3K9Me3ChIP3_readsCounts_GB_SCALED.RDS")

# ORDER QUANTIF decreasingly
Q_eGFPdsRNA_H3_ip_Rep1_GB = Q_eGFPdsRNA_H3_ip_Rep1_GB[order(Q_eGFPdsRNA_H3_ip_Rep1_GB, decreasing=T),drop=F]
Q_GFPRNAi_H3K9Me3ChIP3_GB = Q_GFPRNAi_H3K9Me3ChIP3_GB[order(Q_GFPRNAi_H3K9Me3ChIP3_GB, decreasing=T),drop=F]

# ADD in QUANTIFandSCORES_LIST
QUANTIFandSCORES_LIST$Q_eGFPdsRNA_H3_ip_Rep1_GB = Q_eGFPdsRNA_H3_ip_Rep1_GB
QUANTIFandSCORES_LIST$Q_GFPRNAi_H3K9Me3ChIP3_GB = Q_GFPRNAi_H3K9Me3ChIP3_GB

##################################################
## H3K27ac / H3K4me1 / H3K4me3 (GSE85191)
##################################################
# GNref = readRDS(paste0(workdir, "PROJET_K27K9K36/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

# LOAD QUANTIF
# Q_H3K27ac_rep1_GB = readRDS("/home/depierre/Bureau/work_cperrois/ALIGN_GSE85191/QUANITF/GB/H3K27ac_rep1_readsCounts_GB_SCALED.RDS")
# Q_H3K4me1_rep1_GB = readRDS("/home/depierre/Bureau/work_cperrois/ALIGN_GSE85191/QUANITF/GB/H3K4me1_rep1_readsCounts_GB_SCALED.RDS")
# Q_H3K4me3_rep1_GB = readRDS("/home/depierre/Bureau/work_cperrois/ALIGN_GSE85191/QUANITF/GB/H3K4me3_rep1_readsCounts_GB_SCALED.RDS")

# ORDER QUANTIF decreasingly
Q_H3K27ac_rep1_GB = Q_H3K27ac_rep1_GB[order(Q_H3K27ac_rep1_GB, decreasing=T),drop=F]
Q_H3K4me1_rep1_GB = Q_H3K4me1_rep1_GB[order(Q_H3K4me1_rep1_GB, decreasing=T),drop=F]
Q_H3K4me3_rep1_GB = Q_H3K4me3_rep1_GB[order(Q_H3K4me3_rep1_GB, decreasing=T),drop=F]

# ADD in QUANTIFandSCORES_LIST
QUANTIFandSCORES_LIST$Q_H3K27ac_rep1_GB = Q_H3K27ac_rep1_GB
QUANTIFandSCORES_LIST$Q_H3K4me1_rep1_GB = Q_H3K4me1_rep1_GB
QUANTIFandSCORES_LIST$Q_H3K4me3_rep1_GB = Q_H3K4me3_rep1_GB


#########################################################################################################################
## QUANTIF H3K27me3 GB mutant BEAF (used in FIGURE_6_B_C_HEATMAP_FISHER)
#########################################################################################################################

# LOAD QUANTIF
# Q_H3K27me3_CTL_exp2_GSM3733917 = readRDS("/path/to/QUANITF/GB/H3K27me3_CTL_exp2_GSM3733917_readsCounts_GB.RDS")
# Q_H3K27me3_muBeaf_exp2_GSM3733919 = readRDS("/path/to/QUANITF/GB/H3K27me3_muBeaf_exp2_GSM3733919_readsCounts_GB.RDS")

## COMPUTE ZSCORE
ZSCORE_H3K27me3_muBeaf_exp2 = computeZscore(Q_H3K27me3_muBeaf_exp2_GSM3733919, Q_H3K27me3_CTL_exp2_GSM3733917)

# ORDER ZSCORE decreasingly
ZSCORE_H3K27me3_muBeaf_exp2 = ZSCORE_H3K27me3_muBeaf_exp2[order(ZSCORE_H3K27me3_muBeaf_exp2, decreasing=T),drop=F]

QUANTIFandSCORES_LIST$ZSCORE_H3K27me3_muBeaf_exp2 = ZSCORE_H3K27me3_muBeaf_exp2


#########################################################################################################################
## GENES CP190 BOUND (from GSE52887) on promoter 
#########################################################################################################################
gene_dm6_gr_PROM = promoters(gene_dm6_gr, upstream=500, downstream=250) # Donc je reprend plus large

# LOAD CP190 PEAK CALLING BED 
CP190_1_WT_S2_inputNorm_peaks_1kb = readRDS("/path/to/PEAKCALLING_MACS2/CP190_1_WT_S2_inputNorm_peaks_1kb.RDS")
CP190_2_WT_S2_inputNorm_peaks_1kb = readRDS("/path/to/PEAKCALLING_MACS2/CP190_2_WT_S2_inputNorm_peaks_1kb.RDS")


FindOverlaps = function(GN_GR, peak_GR){
	# seqlevelsStyle(GN_GR) = "UCSC"
	ol1 = findOverlaps(GN_GR, peak_GR, , minoverlap=100,  ignore.strand=TRUE)
	ol1_GN = unique(ol1@from)
	return(GN_GR[ol1_GN])
}

genes_PROM_CP190_1_WT_S2_inputNorm = FindOverlaps(gene_dm6_gr_PROM, CP190_1_WT_S2_inputNorm_peaks_1kb)
genes_PROM_CP190_2_WT_S2_inputNorm = FindOverlaps(gene_dm6_gr_PROM, CP190_2_WT_S2_inputNorm_peaks_1kb)

GENES_CP190 = unique(c(names(genes_PROM_CP190_1_WT_S2_inputNorm), names(genes_PROM_CP190_2_WT_S2_inputNorm)))

QUANTIFandSCORES_LIST$GENES_CP190 = GENES_CP190



#########################################################################################################################
## COMPUTE DISTANCE BETWEEN GENES TSS AND CLOSEST H3K27me3 DOMAINS BORDER
#########################################################################################################################

# LOAD H3K27me3 domains border defined in DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/ folder
K27DOM_BORDER_sup1500 = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed"))
K27DOM_BORDER_sup1500_centered1_oriented <- resize(K27DOM_BORDER_sup1500,1,"center")
K27DOM_BORDER_sup1500_centered1_oriented
# Get TSS coordinates
gene_dm6_gr_TSS = resize(gene_dm6_gr, 1, "start")

# Get distance to nearest for each TSS <-> K27DOM_BORDER
dist_gene2border_TABLE = distanceToNearest(gene_dm6_gr_TSS, K27DOM_BORDER_sup1500_centered1_oriented, ignore.strand=TRUE)

dist_gene2border = dist_gene2border_TABLE@elementMetadata$distance
names(dist_gene2border) = names(gene_dm6_gr_TSS[dist_gene2border_TABLE@from])


QUANTIFandSCORES_LIST$dist_gene2border = dist_gene2border


#########################################################################################################################
# RNASEQ - DIFF EXPR GENES
#########################################################################################################################

LOGFC_DiffExpr_HYPBKD_RPGC_pval005 = readRDS(paste0("/path/to/DATA/DIFF_EXPR_ANALYSIS/DIFF_EXPR_RPGC/LOGFC_DiffExpr_HYPBKD_RPGC_pval005.RDS"))
LOGFC_DiffExpr_MES4KD_RPGC_pval005 = readRDS(paste0("/path/to/DATA/DIFF_EXPR_ANALYSIS/DIFF_EXPR_RPGC/LOGFC_DiffExpr_MES4KD_RPGC_pval005.RDS"))
LOGFC_HYPBKD_RPGC_pval005 = LOGFC_DiffExpr_HYPBKD_RPGC_pval005[,1]
LOGFC_MES4KD_RPGC_pval005 = LOGFC_DiffExpr_MES4KD_RPGC_pval005[,1]
names(LOGFC_HYPBKD_RPGC_pval005) = rownames(LOGFC_DiffExpr_HYPBKD_RPGC_pval005)
names(LOGFC_MES4KD_RPGC_pval005) = rownames(LOGFC_DiffExpr_MES4KD_RPGC_pval005)

LOGFC_HYPBKD_RPGC_pval005 = LOGFC_HYPBKD_RPGC_pval005[order(LOGFC_HYPBKD_RPGC_pval005, decreasing=T),drop=F]
LOGFC_MES4KD_RPGC_pval005 = LOGFC_MES4KD_RPGC_pval005[order(LOGFC_MES4KD_RPGC_pval005, decreasing=T),drop=F]

QUANTIFandSCORES_LIST$LOGFC_HYPBKD_RPGC_pval005 = LOGFC_HYPBKD_RPGC_pval005
QUANTIFandSCORES_LIST$LOGFC_MES4KD_RPGC_pval005 = LOGFC_MES4KD_RPGC_pval005



## ACTIVE GENES FROM RNAseq

ACTIVE_GENES = readRDS(paste0("/path/to/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))
QUANTIFandSCORES_LIST$ACTIVE_GENES = ACTIVE_GENES

#########################################################################################################################
# SAVE QUANTIFandSCORES_LIST
#########################################################################################################################


saveRDS(QUANTIFandSCORES_LIST, paste0(workdir, "DATA_PROCESSING/QUANTIFandSCORES_LIST.RDS"))



#end

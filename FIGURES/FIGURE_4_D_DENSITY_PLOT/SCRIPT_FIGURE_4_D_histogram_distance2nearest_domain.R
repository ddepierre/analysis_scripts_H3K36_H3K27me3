######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# SCRIPT PLOT HSITOGRAM FIGURE 2 C
# Histogram showing distribution density of bins with H3K27me3 increase from H3K27me3 domains border

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

### LOAD histogram funtion
source(paste0(workdir, "FIGURES/R_PLOT_FUNCTION/HISTOGRAM_function.R"))


compute_DIST_feat2BorderOriented.v220503 = function(FEAT1, REFfrom){
  FEAT1 <- resize(FEAT1,1,"center")
  DIST = distanceToNearest(FEAT1, REFfrom, ignore.strand=TRUE)
  dist2border_oriented = unlist(lapply(1:length(DIST@from), function(i){
      # print(i)
      if(as.vector(strand(REFfrom[DIST@to[i]])) %in% "-"){
        start(REFfrom[DIST@to[i]])-start(FEAT1[DIST@from[i]])
      }else{
        start(FEAT1[DIST@from[i]])-start(REFfrom[DIST@to[i]])
      }
  }))
  names(dist2border_oriented) = names(FEAT1[DIST@from])
  return(dist2border_oriented)
}


######################################################################################################################################################
### LOAD DATA
######################################################################################################################################################
### LOAD QUANTIF AND SCORES
QUANTIFandSCORES_LIST = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIFandSCORES_LIST.RDS"))


### LOAD BINS WITH H3K27me3 CALLED WITH NORMR in DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/
DIFFCALLING_K27H_K27C = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_K27H_K27C_bin2000_fdr0001.bed"))
DIFFCALLING_H27_2_C27_1 = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_H27_2_C27_1_bin2000_fdr0001.bed"))
DIFFCALLING_M27_2_C27_1 = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_M27_2_C27_1_bin2000_fdr0001.bed"))
DIFFCALLING_K27M_K27C = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_K27M_K27C_bin2000_fdr0001.bed"))

names(DIFFCALLING_K27H_K27C) = paste0(seqnames(DIFFCALLING_K27H_K27C), "_", start(DIFFCALLING_K27H_K27C))
names(DIFFCALLING_H27_2_C27_1) = paste0(seqnames(DIFFCALLING_H27_2_C27_1), "_", start(DIFFCALLING_H27_2_C27_1))
names(DIFFCALLING_M27_2_C27_1) = paste0(seqnames(DIFFCALLING_M27_2_C27_1), "_", start(DIFFCALLING_M27_2_C27_1))
names(DIFFCALLING_K27M_K27C) = paste0(seqnames(DIFFCALLING_K27M_K27C), "_", start(DIFFCALLING_K27M_K27C))

### LOAD RANDOM 2KB BINS as control
Dmelanogaster_UCSC_dm6_GR_bin2000 = sample(readRDS(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/Dmelanogaster_UCSC_dm6_GR_bin2000.RDS")), 5000)

seqlevelsStyle(DIFFCALLING_K27H_K27C) ="Ensembl"
seqlevelsStyle(DIFFCALLING_H27_2_C27_1) ="Ensembl"
seqlevelsStyle(DIFFCALLING_M27_2_C27_1) ="Ensembl"
seqlevelsStyle(DIFFCALLING_K27M_K27C) ="Ensembl"
seqlevelsStyle(Dmelanogaster_UCSC_dm6_GR_bin2000) ="Ensembl"


### LOAD H3K27me3 DOMAINS BORDER DONE in DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/
H3K27me3_CTRL_DOMAINS_borders = readRDS(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINS_borders_gapsup1500bp_sup1500bp_filledInputGaps.RDS"))

gene_dm6_gr <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
names(gene_dm6_gr) = paste0(names(gene_dm6_gr), ".1")
seqlevelsStyle(gene_dm6_gr) ="Ensembl"

CLUSTERING_DEA_genes = readRDS(paste0(workdir, "FIGURES/FIGURE_3_B_CLUSTERING_HEATMAP/CLUSTERING_DEA_genes.RDS"))


DEA_Cluster_7G_red = list(
											D1D4 = c(CLUSTERING_DEA_genes$D1, CLUSTERING_DEA_genes$D4),
											D2D3 = c(CLUSTERING_DEA_genes$D2, CLUSTERING_DEA_genes$D3))




######################################################################################################################################################
### Compute For each bin with H3K27me3 increase the distance to nearest domains border
######################################################################################################################################################

OVLP_DIFFCALLING_K27H_K27C_DEA_D1D4 = findOverlaps(resize(DIFFCALLING_K27H_K27C, 2000, "center"), gene_dm6_gr[DEA_Cluster_7G_red$D1D4], ignore.strand=T, type="any", select="all", maxgap=2000)
DIFFCALLING_K27H_K27C_DEA_D1D4 = DIFFCALLING_K27H_K27C[unique(OVLP_DIFFCALLING_K27H_K27C_DEA_D1D4@from)]

OVLP_DIFFCALLING_K27M_K27C_DEA_D2D3 = findOverlaps(resize(DIFFCALLING_K27M_K27C, 2000, "center"), gene_dm6_gr[DEA_Cluster_7G_red$D2D3], ignore.strand=T, type="any", select="all", maxgap=2000)
DIFFCALLING_K27M_K27C_DEA_D2D3 = DIFFCALLING_K27M_K27C[unique(OVLP_DIFFCALLING_K27M_K27C_DEA_D2D3@from)]


DISTfromBORDER_DIFFCALLING_K27H_K27C_DEA_D1D4 = compute_DIST_feat2BorderOriented.v220503(FEAT1 = DIFFCALLING_K27H_K27C_DEA_D1D4, REFfrom = H3K27me3_CTRL_DOMAINS_borders)
DISTfromBORDER_DIFFCALLING_K27M_K27C_DEA_D2D3 = compute_DIST_feat2BorderOriented.v220503(FEAT1 = DIFFCALLING_K27M_K27C_DEA_D2D3, REFfrom = H3K27me3_CTRL_DOMAINS_borders)
DISTfromBORDER_Dmelanogaster_UCSC_dm6_GR_bin2000 = compute_DIST_feat2BorderOriented.v220503(FEAT1 = Dmelanogaster_UCSC_dm6_GR_bin2000, REFfrom = H3K27me3_CTRL_DOMAINS_borders)
names(DISTfromBORDER_Dmelanogaster_UCSC_dm6_GR_bin2000) = paste0("random_",names(DISTfromBORDER_Dmelanogaster_UCSC_dm6_GR_bin2000))


######################################################################################################################################################
### plot histogram density of distances
######################################################################################################################################################


VecNames = names(DISTfromBORDER_DIFFCALLING_K27H_K27C_DEA_D1D4)
dist_bin2border = c(DISTfromBORDER_DIFFCALLING_K27H_K27C_DEA_D1D4, DISTfromBORDER_Dmelanogaster_UCSC_dm6_GR_bin2000)
plotHisto_listGN(vec = dist_bin2border, nameVec = "dist_bin2border", vec_GN = VecNames, nameVec_GN = "DISTfromBORDER_DIFFCALLING_K27H_K27C_DEA_D1D4", xlim = c(-20000, 20000), ylim = c(0,0.0002), binW = 2000, alphaTransp = 0.6,
                fillingcol = c("#3d3d3d", "#820002"),outDir = paste0(workdir, "FIGURES/FIGURE_4_D_DENSITY_PLOT/"), info = "FIGURE_4_D")


VecNames = names(DISTfromBORDER_DIFFCALLING_K27M_K27C_DEA_D2D3)
dist_bin2border = c(DISTfromBORDER_DIFFCALLING_K27M_K27C_DEA_D2D3, DISTfromBORDER_Dmelanogaster_UCSC_dm6_GR_bin2000)
plotHisto_listGN(vec = dist_bin2border, nameVec = "dist_bin2border", vec_GN = VecNames, nameVec_GN = "DISTfromBORDER_DIFFCALLING_K27M_K27C_DEA_D2D3", xlim = c(-20000, 20000), ylim = c(0,0.0002), binW = 2000, alphaTransp = 0.6,
                fillingcol = c("#3d3d3d", "#820002"),outDir = paste0(workdir, "FIGURES/FIGURE_4_D_DENSITY_PLOT/"), info = "FIGURE_4_D")


######################################################################################################################################################
### END
######################################################################################################################################################



#end

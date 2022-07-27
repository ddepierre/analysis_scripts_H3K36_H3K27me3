######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# QUANTIFICATION ChIP-seq on Gene Body

######################################################################################################################################################
######################################################################################################################################################
# R version 3.4.2 (2017-09-28) -- "Short Summer"
# Copyright (C) 2017 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

### LIBRARY

######################################################################################################################################################
######################################################################################################################################################
### LOAD DATA
workdir = ""

prof_dir = paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/")
quantif_dir = paste0(workdir, "DATA_PROCESSING/QUANTIF/GENE_BODY/")


# Load data
C27_1_RPGC = readRDS(paste0(prof_dir, "C27_1_trimmed_filt_sort_RPGC_profmat.RDS"))
H27_2_RPGC = readRDS(paste0(prof_dir, "H27_2_trimmed_filt_sort_RPGC_profmat.RDS"))
H3K36me2_RPGC = readRDS(paste0(prof_dir, "H3K36me2_chipseq_filt_sort_RPGC_profmat.RDS"))
H3K36me3_RPGC = readRDS(paste0(prof_dir, "H3K36me3_chipseq_filt_sort_RPGC_profmat.RDS"))
K27C_RPGC = readRDS(paste0(prof_dir, "K27C_trimmed_filt_sort_RPGC_profmat.RDS"))
K27H_RPGC = readRDS(paste0(prof_dir, "K27H_trimmed_filt_sort_RPGC_profmat.RDS"))
K27M_RPGC = readRDS(paste0(prof_dir, "K27M_trimmed_filt_sort_RPGC_profmat.RDS"))
M27_2_RPGC = readRDS(paste0(prof_dir, "M27_2_trimmed_filt_sort_RPGC_profmat.RDS"))



# get reads count on a given window
readSumWindow = function(prof, start=400, end=800){
	Nameprof = deparse(substitute(prof))
	print(Nameprof)
	prof = rowSums(prof[,c(start:end)])
	filenameRDS = paste0(quantif_dir,Nameprof ,  "_readsCounts_GB_SCALED.RDS")
	saveRDS(prof, file=filenameRDS)
}


readSumWindow(C27_1_RPGC)
readSumWindow(H27_2_RPGC)
readSumWindow(H3K36me2_RPGC)
readSumWindow(H3K36me3_RPGC)
readSumWindow(K27C_RPGC)
readSumWindow(K27H_RPGC)
readSumWindow(K27M_RPGC)
readSumWindow(M27_2_RPGC)

Q_C27_1_RPGC = readRDS(paste0(quantif_dir, "C27_1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H27_2_RPGC = readRDS(paste0(quantif_dir, "H27_2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me2_RPGC = readRDS(paste0(quantif_dir, "H3K36me2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me3_RPGC = readRDS(paste0(quantif_dir, "H3K36me3_RPGC_readsCounts_GB_SCALED.RDS"))
Q_K27C_RPGC = readRDS(paste0(quantif_dir, "K27C_RPGC_readsCounts_GB_SCALED.RDS"))
Q_K27H_RPGC = readRDS(paste0(quantif_dir, "K27H_RPGC_readsCounts_GB_SCALED.RDS"))
Q_K27M_RPGC = readRDS(paste0(quantif_dir, "K27M_RPGC_readsCounts_GB_SCALED.RDS"))
Q_M27_2_RPGC = readRDS(paste0(quantif_dir, "M27_2_RPGC_readsCounts_GB_SCALED.RDS"))



#end

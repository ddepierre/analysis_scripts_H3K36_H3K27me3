######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# QUANTIFICATION of CHIPSEQ READS at the H3K27me3_DOMAINS_BORDER

######################################################################################################################################################
######################################################################################################################################################
# R version 3.4.2 (2017-09-28) -- "Short Summer"
# Copyright (C) 2017 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

### LIBRARY

######################################################################################################################################################
######################################################################################################################################################
### LOAD DATA
workdir = "" # on CUVIER10

prof_dir = paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/")
quantif_dir = paste0(workdir, "DATA_PROCESSING/QUANTIF/H3K27me3_DOMAINS_BORDER/")
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
readSumWindow = function(prof, start=450, end=550){
	Nameprof = deparse(substitute(prof))
	print(Nameprof)
	prof = rowSums(prof[,c(start:end)])
	filenameRDS = paste0(quantif_dir,Nameprof ,  "_readsCounts_BORDER_1kb.RDS")
	saveRDS(prof, file=filenameRDS)
}

readSumWindow(H3K36me2_RPGC)
readSumWindow(H3K36me3_RPGC)

Q_H3K36me2_RPGC = readRDS(paste0(quantif_dir, "H3K36me2_RPGC_readsCounts_BORDER_1kb.RDS"))
Q_H3K36me3_RPGC = readRDS(paste0(quantif_dir, "H3K36me3_RPGC_readsCounts_BORDER_1kb.RDS"))


#end

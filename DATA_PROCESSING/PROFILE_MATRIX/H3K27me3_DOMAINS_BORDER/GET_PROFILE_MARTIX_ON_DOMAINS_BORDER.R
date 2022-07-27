######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# GET PROFILE MATRIX ON GENE BODY FOR EACH GENES

######################################################################################################################################################
######################################################################################################################################################


saveRDS(K27DOM_BORDER_sup6000, paste0(workdir, "DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup6000_BORDER_centered_oriented.RDS"))
export.bed(K27DOM_BORDER_sup6000,paste0(workdir,"DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup6000_BORDER_centered_oriented.bed"))






module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/C27_1_trimmed_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/C27_1_trimmed_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/C27_1_trimmed_filt_sort_RPGC.outmatrix
module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/H27_2_trimmed_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H27_2_trimmed_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H27_2_trimmed_filt_sort_RPGC.outmatrix
module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/M27_2_trimmed_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/M27_2_trimmed_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/M27_2_trimmed_filt_sort_RPGC.outmatrix

module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/K27C_trimmed_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27C_trimmed_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27C_trimmed_filt_sort_RPGC.outmatrix
module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/K27H_trimmed_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27H_trimmed_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27H_trimmed_filt_sort_RPGC.outmatrix
module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/K27M_trimmed_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27M_trimmed_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27M_trimmed_filt_sort_RPGC.outmatrix
# module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/Input27_1_trimmed_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/Input27_1_trimmed_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/Input27_1_trimmed_filt_sort_RPGC.outmatrix
# module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/Input_trimmed_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/Input_trimmed_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/Input_trimmed_filt_sort_RPGC.outmatrix
module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/H3K36me2_chipseq_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H3K36me2_chipseq_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H3K36me2_chipseq_filt_sort_RPGC.outmatrix
module load bioinfo/deepTools-3.0.2-python-3.4.3; computeMatrix reference-point -S /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROCESSED_BIGWIG_FILES/H3K36me3_chipseq_filt_sort_RPGC.bw -R /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed -bs 10 --referencePoint center -a 5000 -b 5000 -out /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H3K36me3_chipseq_filt_sort_RPGC.out --outFileNameMatrix /home/ddepierre/work/PROJET_K27K9K36/DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H3K36me3_chipseq_filt_sort_RPGC.outmatrix






########################################################
## 1/ deepTools computeMatrix (BASH)
########################################################
# deepTools/deepTools-2.5.3
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/C27_1_trimmed_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/C27_1_trimmed_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/C27_1_trimmed_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/H27_2_trimmed_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/H27_2_trimmed_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/H27_2_trimmed_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/H3K36me2_chipseq_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/H3K36me2_chipseq_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/H3K36me2_chipseq_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/H3K36me3_chipseq_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/H3K36me3_chipseq_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/H3K36me3_chipseq_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/Input27_1_trimmed_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/Input27_1_trimmed_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/Input27_1_trimmed_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/Input_trimmed_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/Input_trimmed_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/Input_trimmed_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/K27C_trimmed_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/K27C_trimmed_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/K27C_trimmed_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/K27H_trimmed_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/K27H_trimmed_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/K27H_trimmed_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/K27M_trimmed_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/K27M_trimmed_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/K27M_trimmed_filt_sort_RPGC.outmatrix
computeMatrix scale-regions -S DATA_PROCESSING/PROCESSED_BIGWIG_FILES/M27_2_trimmed_filt_sort_RPGC.bw -R DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.bed -bs 10 -m 4000 -a 2000 -b 4000 -out DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/M27_2_trimmed_filt_sort_RPGC.out --outFileNameMatrix DATA_PROCESSING/PROFILE_MATRIX/GENE_BODY/M27_2_trimmed_filt_sort_RPGC.outmatrix

########################################################
## 2/ ADD GENE NAMES ON MATRIX ON SAVE THEM IN R FORMAT
########################################################
# R version 3.4.2 (2017-09-28) -- "Short Summer"
# Copyright (C) 2017 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

### LIBRARY
library("GenomicRanges")
library("S4Vectors")
library("stats4")
library("IRanges")
library("GenomeInfoDb")
library("BiocGenerics")
library("parallel")


######################################################################################################################################################
######################################################################################################################################################
### LOAD DATA
workdir = ""

ref_border_gr = readRDS(paste0(workdir, "DATA_PROCESSING/NOMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.RDS"))

file_outmatrix = as.matrix(read.table(paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27C_trimmed_filt_sort_RPGC.outmatrix"), skip=3))
out_profmat_RDS = paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27C_trimmed_filt_sort_RPGC_profmat")
### ADD GENES ID to as rownames
rownames(file_outmatrix) = paste0(names(ref_border_gr))
saveRDS(file_outmatrix, paste0(out_profmat_RDS,".RDS"))

file_outmatrix = as.matrix(read.table(paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27H_trimmed_filt_sort_RPGC.outmatrix"), skip=3))
out_profmat_RDS = paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27H_trimmed_filt_sort_RPGC_profmat")
### ADD GENES ID to as rownames
rownames(file_outmatrix) = paste0(names(ref_border_gr))
saveRDS(file_outmatrix, paste0(out_profmat_RDS,".RDS"))

file_outmatrix = as.matrix(read.table(paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27M_trimmed_filt_sort_RPGC.outmatrix"), skip=3))
out_profmat_RDS = paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/K27M_trimmed_filt_sort_RPGC_profmat")
### ADD GENES ID to as rownames
rownames(file_outmatrix) = paste0(names(ref_border_gr))
saveRDS(file_outmatrix, paste0(out_profmat_RDS,".RDS"))

file_outmatrix = as.matrix(read.table(paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H3K36me2_chipseq_filt_sort_RPGC.outmatrix"), skip=3))
out_profmat_RDS = paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H3K36me2_chipseq_filt_sort_RPGC_profmat")
### ADD GENES ID to as rownames
rownames(file_outmatrix) = paste0(names(ref_border_gr))
saveRDS(file_outmatrix, paste0(out_profmat_RDS,".RDS"))

file_outmatrix = as.matrix(read.table(paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H3K36me3_chipseq_filt_sort_RPGC.outmatrix"), skip=3))
out_profmat_RDS = paste0(workdir, "DATA_PROCESSING/PROFILE_MATRIX/H3K27me3_DOMAINS_BORDER/H3K36me3_chipseq_filt_sort_RPGC_profmat")
### ADD GENES ID to as rownames
rownames(file_outmatrix) = paste0(names(ref_border_gr))
saveRDS(file_outmatrix, paste0(out_profmat_RDS,".RDS"))



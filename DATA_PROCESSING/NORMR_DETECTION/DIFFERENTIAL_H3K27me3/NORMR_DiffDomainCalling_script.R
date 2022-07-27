######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# DETECTION of H3K27me3 domains with NORMR

######################################################################################################################################################
######################################################################################################################################################
# R version 3.4.2 (2017-09-28) -- "Short Summer"
# Copyright (C) 2017 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

### LIBRARY
library("dplyr")
library("Rsamtools")
library("GenomicRanges")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("normr")
library("BSgenome.Dmelanogaster.UCSC.dm6")
library("rtracklayer")
'%ni%' = Negate('%in%')

######################################################################################################################################################
######################################################################################################################################################
### LOAD DATA
workdir = ""

### BAM PATH
bamdir = paste0(workdir, "DATA_PROCESSING/PROCESSED_BAM_FILES/")

H3K27me3_CTRL_C27_1 <- paste0(bamdir, "C27_1_trimmed_filt_sort.bam")
H3K27me3_HYPB_H27_2 <- paste0(bamdir, "H27_2_trimmed_filt_sort.bam")
H3K27me3_MES4_M27_2 <- paste0(bamdir, "M27_2_trimmed_filt_sort.bam")

H3K27me3_CTRL_K27C <- paste0(bamdir, "K27C_trimmed_filt_sort.bam")
H3K27me3_HYPB_K27H <- paste0(bamdir, "K27H_trimmed_filt_sort.bam")
H3K27me3_MES4_K27M <- paste0(bamdir, "K27M_trimmed_filt_sort.bam")


## LOAD Chrom infos
ChromSize_dm6 <- getChromInfoFromUCSC("dm6")
ChromSize_dm6 <- ChromSize_dm6[1:7,]
ChromSize_dm6$chrom = c("2L", "2R", "3L", "3R", "4",  "X",  "Y")


############################################################################################################################################################################################################################################################
###########################################################                 DIFFERENTIAL CALLING                              ##############################################################################################################################
############################################################################################################################################################################################################################################################
####  ENRICHMENT CALLING -> bin2000
countConfiguration <- countConfigSingleEnd(binsize = 2000 ,mapq = 30, shift = 5)

################## HYPB treatment = H3K27me3_HYPB_K27H, control = H3K27me3_CTRL_K27C
enrich <- enrichR(treatment = H3K27me3_HYPB_K27H, control = H3K27me3_CTRL_K27C, genome = gr,  countConfig = countConfiguration,iterations = 10, procs = 1, verbose = TRUE)
exportR(enrich, filename = paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_K27H_K27C_bin2000.bw"))
GR_enrich_K27_fdr = getRanges(enrich, fdr = 0.001)
seqlevelsStyle(GR_enrich_K27_fdr) = "UCSC"
names(GR_enrich_K27_fdr) = paste0(as.vector(seqnames(GR_enrich_K27_fdr)), "_", start(GR_enrich_K27_fdr))
GR_enrich_K27_fdr # 1008
rtracklayer::export.bed(GR_enrich_K27_fdr,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_K27H_K27C_bin2000_fdr0001.bed"))
saveRDS(GR_enrich_K27_fdr, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_K27H_K27C_bin2000_fdr0001.RDS"))

################## HYPB treatment = H3K27me3_HYPB_H27_2, control = H3K27me3_CTRL_C27_1
enrich <- enrichR(treatment = H3K27me3_HYPB_H27_2, control = H3K27me3_CTRL_C27_1, genome = gr,  countConfig = countConfiguration,iterations = 10, procs = 1, verbose = TRUE)
exportR(enrich, filename = paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_H27_2_C27_1_bin2000.bw"))
GR_enrich_K27_fdr = getRanges(enrich, fdr = 0.001)
seqlevelsStyle(GR_enrich_K27_fdr) = "UCSC"
names(GR_enrich_K27_fdr) = paste0(as.vector(seqnames(GR_enrich_K27_fdr)), "_", start(GR_enrich_K27_fdr))
GR_enrich_K27_fdr  # 179
rtracklayer::export.bed(GR_enrich_K27_fdr,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_H27_2_C27_1_bin2000_fdr0001.bed"))
saveRDS(GR_enrich_K27_fdr, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_H27_2_C27_1_bin2000_fdr0001.RDS"))

######################################## MES4
################## MES4 treatment = H3K27me3_MES4_M27_2, control = H3K27me3_CTRL_C27_1
enrich <- enrichR(treatment = H3K27me3_MES4_M27_2, control = H3K27me3_CTRL_C27_1, genome = gr,  countConfig = countConfiguration,iterations = 10, procs = 1, verbose = TRUE)
exportR(enrich, filename = paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_M27_2_C27_1_bin2000.bw"))
GR_enrich_K27_fdr = getRanges(enrich, fdr = 0.001)
seqlevelsStyle(GR_enrich_K27_fdr) = "UCSC"
names(GR_enrich_K27_fdr) = paste0(as.vector(seqnames(GR_enrich_K27_fdr)), "_", start(GR_enrich_K27_fdr))
GR_enrich_K27_fdr # 1630
rtracklayer::export.bed(GR_enrich_K27_fdr,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_M27_2_C27_1_bin2000_fdr0001.bed"))
saveRDS(GR_enrich_K27_fdr, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_M27_2_C27_1_bin2000_fdr0001.RDS"))

################## MES4 treatment = H3K27me3_MES4_K27M, control = H3K27me3_CTRL_K27C
enrich <- enrichR(treatment = H3K27me3_MES4_K27M, control = H3K27me3_CTRL_K27C, genome = gr,  countConfig = countConfiguration,iterations = 10, procs = 1, verbose = TRUE)
exportR(enrich, filename = paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_K27M_K27C_bin2000.bw"))
GR_enrich_K27_fdr = getRanges(enrich, fdr = 0.001)
seqlevelsStyle(GR_enrich_K27_fdr) = "UCSC"
names(GR_enrich_K27_fdr) = paste0(as.vector(seqnames(GR_enrich_K27_fdr)), "_", start(GR_enrich_K27_fdr))
GR_enrich_K27_fdr # 5288
rtracklayer::export.bed(GR_enrich_K27_fdr,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_K27M_K27C_bin2000_fdr0001.bed"))
saveRDS(GR_enrich_K27_fdr, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/DIFFCALLING_K27M_K27C_bin2000_fdr0001.RDS"))

################################################################################################################################################################
######################################## GET RANDOM BIN FROM GENOME at the same size then DIFFERENTIAL DETECTED BINS ABOVE

Dmelanogaster_UCSC_dm6_GR = GRanges(seqinfo(BSgenome.Dmelanogaster.UCSC.dm6))
Dmelanogaster_UCSC_dm6_GR <- Dmelanogaster_UCSC_dm6_GR[seqnames(Dmelanogaster_UCSC_dm6_GR) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4",  "chrX",  "chrY")]
Dmelanogaster_UCSC_dm6_GR_bin2000 = tile(Dmelanogaster_UCSC_dm6_GR, width=2000)
Dmelanogaster_UCSC_dm6_GR_bin2000 = unlist(Dmelanogaster_UCSC_dm6_GR_bin2000)
names(Dmelanogaster_UCSC_dm6_GR_bin2000) = paste0(as.vector(seqnames(Dmelanogaster_UCSC_dm6_GR_bin2000)), "_", start(Dmelanogaster_UCSC_dm6_GR_bin2000))
saveRDS(Dmelanogaster_UCSC_dm6_GR_bin2000, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DIFFERENTIAL_H3K27me3/Dmelanogaster_UCSC_dm6_GR_bin2000.RDS"))


#end

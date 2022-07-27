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
### LOAD DATA
######################################################################################################################################################
 
workdir = "/home/depierre/Bureau/analysis_scripts_H3K36_H3K27me3/" 

### BAM PATH
bamdir = paste0(workdir, "DATA_PROCESSING/PROCESSED_BAM_FILES/")

H3K27me3_INPUT_Input27_1 <- paste0(bamdir, "Input27_1_trimmed_filt_sort.bam")
H3K27me3_CTRL_C27_1 <- paste0(bamdir, "C27_1_trimmed_filt_sort.bam")


ChromSize_dm6 <- getChromInfoFromUCSC("dm6")
ChromSize_dm6 <- ChromSize_dm6[1:7,]
ChromSize_dm6$chrom = c("2L", "2R", "3L", "3R", "4",  "X",  "Y")

###########################
## BIN de 200
###########################
countConfiguration <- countConfigSingleEnd(binsize = 200 ,mapq = 30, shift = 20)

enrich <- enrichR(treatment = H3K27me3_CTRL_C27_1, control = H3K27me3_INPUT_Input27_1, genome = ChromSize_dm6,  countConfig = countConfiguration,iterations = 10, procs = 1, verbose = TRUE)
exportR(enrich, filename = paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_norminput_bin200.bw"))
GR_enrich_K27_fdr = getRanges(enrich, fdr = 0.0001)
seqlevelsStyle(GR_enrich_K27_fdr) = "UCSC"

GR_enrich_K27_fdr_reduced = reduce(GR_enrich_K27_fdr, drop.empty.ranges=FALSE, min.gapwidth=1500,with.revmap=FALSE, with.inframe.attrib=FALSE)
GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp = GR_enrich_K27_fdr_reduced[width(GR_enrich_K27_fdr_reduced)>1500]

names(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp) = paste0(as.vector(seqnames(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)), "_", start(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp))

rtracklayer::export.bed(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSoverInput_bin200_fdr00001_gapsup1500b_sup1500b.bed"))
saveRDS(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSoverInput_bin200_fdr00001_gapsup1500b_sup1500b.RDS"))







#  create a bed and GR from OUT domain for PROFILE_MATRIX_BORDER

GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp = rtracklayer::import.bed(paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSoverInput_bin200_fdr00001_gapsup1500b_sup1500b.bed"))

K27DOM_BORDER1 = GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp
start(K27DOM_BORDER1) = end(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)-5000
end(K27DOM_BORDER1) = end(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)+5000
strand(K27DOM_BORDER1) = "+"
K27DOM_BORDER2 = GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp
start(K27DOM_BORDER2) = start(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)-5000
end(K27DOM_BORDER2) = start(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)+5000
strand(K27DOM_BORDER2) = "-"
K27DOM_BORDER_sup1500 = c(K27DOM_BORDER1,K27DOM_BORDER2)
names(K27DOM_BORDER_sup1500) = paste0(as.vector(seqnames(K27DOM_BORDER_sup1500)), "_", start(K27DOM_BORDER_sup1500))
K27DOM_BORDER_sup1500$name = names(K27DOM_BORDER_sup1500)

# K27DOM_BORDER_sup1500
# GRanges object with 7770 ranges and 1 metadata column:


saveRDS(K27DOM_BORDER_sup1500, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.RDS"))
rtracklayer::export.bed(K27DOM_BORDER_sup1500,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSsup1500_BORDER_centered_oriented.bed"))




#######################################################################################################################################
# In order to refine H3K27me3 detected domains, we are oing to filter them by null region in input, to avoid considering end of a domain as a border, when it is actually a none mappable region or a null signal region.
#######################################################################################################################################

# 1/ from ChIP-seq input, get none mappable/no signal region, defined as region with score = 0, and with length > 1500bp 
	# chromosomes.chr <- c("2L", "2R", "3L", "3R", "4",  "X",  "Y")
    # Input_trimmed_filt_sort_RPGC.bw = rtracklayer::import.bw(paste0(workdir, "DATA/BIGWIG/Input_trimmed_filt_sort_RPGC.bw"))
    # Input_trimmed_filt_sort_RPGC.bw = Input_trimmed_filt_sort_RPGC.bw[seqnames(Input_trimmed_filt_sort_RPGC.bw) %in% chromosomes.chr]
    # Input_trimmed_filt_sort_RPGC.bw = reduce(Input_trimmed_filt_sort_RPGC.bw[Input_trimmed_filt_sort_RPGC.bw$score == 0])
    # Input_trimmed_filt_sort_RPGC.bw = Input_trimmed_filt_sort_RPGC.bw[width(Input_trimmed_filt_sort_RPGC.bw)>1500]
    # Input27_1_filt_sort_RPGC.bw = rtracklayer::import.bw(paste0(workdir, "DATA/BIGWIG/Input27_1_filt_sort_RPGC.bw"))
    # Input27_1_filt_sort_RPGC.bw = Input27_1_filt_sort_RPGC.bw[seqnames(Input27_1_filt_sort_RPGC.bw) %in% chromosomes.chr]
    # Input27_1_filt_sort_RPGC.bw = reduce(Input27_1_filt_sort_RPGC.bw[Input27_1_filt_sort_RPGC.bw$score == 0])
    # Input27_1_filt_sort_RPGC.bw = Input27_1_filt_sort_RPGC.bw[width(Input27_1_filt_sort_RPGC.bw)>1500]
	# none_mappable_regions_from_input_sup1500 = reduce(c(Input_trimmed_filt_sort_RPGC.bw, Input27_1_filt_sort_RPGC.bw))
    # rtracklayer::export.bed(none_mappable_regions_from_input_sup1500, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/none_mappable_regions_from_input_sup1500_sup1500.bed"))

	nullSignal_region = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/none_mappable_regions_from_input_sup1500.bed"))

# 2/ from nullSignal_region, filter gaps in H3K27me3 domains 
	GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINSoverInput_bin200_fdr00001_gapsup1500b_sup1500b.bed"))
	GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp
	seqlevelsStyle(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp) = "ensembl"

	gaps_K27_dom =gaps(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)
	gaps_K27_dom = gaps_K27_dom[seqnames(gaps_K27_dom) %in% chromosomes.chr]
	gaps_K27_dom = gaps_K27_dom[strand(gaps_K27_dom) %in% "*"]
	gaps_K27_dom_ovlp_nullSignal_region = subsetByOverlaps(gaps_K27_dom, nullSignal_region)

	gaps_K27_dom_ovlp_nullSignal_region = findOverlaps(gaps_K27_dom, nullSignal_region)

    OVLP_gaps_nonMap.gr <- pintersect(gaps_K27_dom[queryHits(gaps_K27_dom_ovlp_nullSignal_region)], nullSignal_region[subjectHits(gaps_K27_dom_ovlp_nullSignal_region)])
    gaps_K27_dom_ovlp_nullSignal_region.gr = gaps_K27_dom[queryHits(gaps_K27_dom_ovlp_nullSignal_region)]
    percentOverlap <- width(OVLP_gaps_nonMap.gr) / width(gaps_K27_dom[queryHits(gaps_K27_dom_ovlp_nullSignal_region)])
	gaps_K27_dom_ovlp_nullSignal_region.gr$percentOverlap = percentOverlap

	gaps_K27_dom_ovlp_nullSignal_region_sup65.gr = gaps_K27_dom_ovlp_nullSignal_region.gr[gaps_K27_dom_ovlp_nullSignal_region.gr$percentOverlap > 0.65]
	
	gaps_K27_dom_ovlp_nullSignal_region_sup65.gr$percentOverlap = NULL
	GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp$name = NULL
	GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp$score = NULL
	GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp_filledInputGaps = reduce(c(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp, gaps_K27_dom_ovlp_nullSignal_region_sup65.gr))
	rtracklayer::export.bed(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp_filledInputGaps, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp_filledInputGaps.bed"))

	# save oriented borders
	H3K27me3_CTRL_DOMAINS_filledGaps = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp_filledInputGaps.bed"))
	seqlevelsStyle(H3K27me3_CTRL_DOMAINS_filledGaps) = 'Ensembl'
    chromosomes.chr <- c('2L','2R','3R','3L','X')
  	
    H3K27me3_CTRL_DOMAINS_filledGaps = H3K27me3_CTRL_DOMAINS_filledGaps[seqnames(H3K27me3_CTRL_DOMAINS_filledGaps) %in% chromosomes.chr]
    H3K27me3_CTRL_DOMAINS_filledGaps = sort(H3K27me3_CTRL_DOMAINS_filledGaps)
  	
    H3K27me3_CTRL_DOMAINS_start = resize(H3K27me3_CTRL_DOMAINS_filledGaps, fix="start", 1)
    strand(H3K27me3_CTRL_DOMAINS_start) = "-"
    H3K27me3_CTRL_DOMAINS_end = resize(H3K27me3_CTRL_DOMAINS_filledGaps, fix="end", 1)
    strand(H3K27me3_CTRL_DOMAINS_end) = "+"
    H3K27me3_CTRL_DOMAINS_borders = c(H3K27me3_CTRL_DOMAINS_start, H3K27me3_CTRL_DOMAINS_end)

  	names(H3K27me3_CTRL_DOMAINS_borders) = paste0("BorDOM_",seqnames(H3K27me3_CTRL_DOMAINS_borders), "_", start(H3K27me3_CTRL_DOMAINS_borders))
  	
    rtracklayer::export.bed(H3K27me3_CTRL_DOMAINS_borders, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINS_borders_gapsup1500bp_sup1500bp_filledInputGaps.bed"))
    saveRDS(H3K27me3_CTRL_DOMAINS_borders, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINS_borders_gapsup1500bp_sup1500bp_filledInputGaps.RDS"))



####################################################################################################################################################################
## LABEL GENES ACCORDING TO THEIR POSITION ON H3K27me3 domains
####################################################################################################################################################################

# create a granges list with 4 classes : H3K27me3_INDOM, H3K27me3_OUTDOM, H3K27me3_INBORDER, H3K27me3_OUTBORDER
H3K27me3_OUTDOM_tmp = gaps(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)
H3K27me3_OUTDOM_tmp = H3K27me3_OUTDOM_tmp[strand(H3K27me3_OUTDOM_tmp) %in% "*"]
H3K27me3_OUTBORDER_tmp1 = H3K27me3_OUTDOM_tmp[width(H3K27me3_OUTDOM_tmp)<5001]
H3K27me3_OUTDOM_tmp = H3K27me3_OUTDOM_tmp[width(H3K27me3_OUTDOM_tmp)>5000]
H3K27me3_OUTDOM = H3K27me3_OUTDOM_tmp
start(H3K27me3_OUTDOM)=start(H3K27me3_OUTDOM)+2500
end(H3K27me3_OUTDOM)=end(H3K27me3_OUTDOM)-2500
H3K27me3_OUTDOM = H3K27me3_OUTDOM[width(H3K27me3_OUTDOM)>500] #2308
names(H3K27me3_OUTDOM) = paste0(as.vector(seqnames(H3K27me3_OUTDOM)), "_", start(H3K27me3_OUTDOM))
# rtracklayer::export.bed(H3K27me3_OUTDOM,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_OUTDOM_official2.bed"))

H3K27me3_OUTBORDER_tmp2 = H3K27me3_OUTDOM_tmp
end(H3K27me3_OUTBORDER_tmp2) = start(H3K27me3_OUTBORDER_tmp2)+2499
H3K27me3_OUTBORDER_tmp3 = H3K27me3_OUTDOM_tmp
start(H3K27me3_OUTBORDER_tmp3) = end(H3K27me3_OUTBORDER_tmp3)-2499
H3K27me3_OUTBORDER = c(H3K27me3_OUTBORDER_tmp1, H3K27me3_OUTBORDER_tmp2, H3K27me3_OUTBORDER_tmp3) #6301
names(H3K27me3_OUTBORDER) = paste0(as.vector(seqnames(H3K27me3_OUTBORDER)), "_", start(H3K27me3_OUTBORDER))
# rtracklayer::export.bed(H3K27me3_OUTBORDER,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_OUTBORDER_official2.bed"))

H3K27me3_INDOM_tmp = GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp
H3K27me3_INBORDER_tmp1 = GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp[width(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)<5001]
H3K27me3_INDOM = GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp[width(GR_enrich_K27_fdr_reduced_gapsup1500bp_sup1500bp)>5000]
start(H3K27me3_INDOM)=start(H3K27me3_INDOM)+2500
end(H3K27me3_INDOM)=end(H3K27me3_INDOM)-2500
H3K27me3_INDOM = H3K27me3_INDOM[width(H3K27me3_INDOM)>200] #2577
names(H3K27me3_INDOM) = paste0(as.vector(seqnames(H3K27me3_INDOM)), "_", start(H3K27me3_INDOM))
# rtracklayer::export.bed(H3K27me3_INDOM,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_INDOM_official2.bed"))

H3K27me3_INBORDER_tmp2 = H3K27me3_INDOM_tmp
end(H3K27me3_INBORDER_tmp2) = start(H3K27me3_INBORDER_tmp2)+2499
H3K27me3_INBORDER_tmp3 = H3K27me3_INDOM_tmp
start(H3K27me3_INBORDER_tmp3) = end(H3K27me3_INBORDER_tmp3)-2499
H3K27me3_INBORDER = c(H3K27me3_INBORDER_tmp1, H3K27me3_INBORDER_tmp2, H3K27me3_INBORDER_tmp3) #6301
names(H3K27me3_INBORDER) = paste0(as.vector(seqnames(H3K27me3_INBORDER)), "_", start(H3K27me3_INBORDER))
# rtracklayer::export.bed(H3K27me3_INBORDER,paste0(workdir,"DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_INBORDER_official2.bed"))

H3K27me3_DOMAIN_GR = list(H3K27me3_OUTDOM, H3K27me3_OUTBORDER, H3K27me3_INDOM, H3K27me3_INBORDER)
names(H3K27me3_DOMAIN_GR) = c("H3K27me3_OUTDOM", "H3K27me3_OUTBORDER", "H3K27me3_INDOM", "H3K27me3_INBORDER")

saveRDS(H3K27me3_DOMAIN_GR, paste0(workdir, "PROJET_K27K9K36/DATA/DOMAIN_NORMR/H3K27me3_DOMAINS_GRList.RDS"))
lapply(H3K27me3_DOMAIN_GR, length)

# lapply(H3K27me3_DOMAIN_GR, length)
# $H3K27me3_OUTDOM
# [1] 2308

# $H3K27me3_OUTBORDER
# [1] 6301

# $H3K27me3_INDOM
# [1] 2577

# $H3K27me3_INBORDER
# [1] 9025



# get GENES in different domains/border
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
gene_dm6_gr <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
names(gene_dm6_gr) = paste0(names(gene_dm6_gr),".1")
gene_dm6_gr <- gene_dm6_gr[seqnames(gene_dm6_gr) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4",  "chrX",  "chrY")]

myfindOverlaps = function(genes_GR, ref_GR){
	# seqlevelsStyle(peak_GR) = "UCSC"
	ol1 = findOverlaps(genes_GR, ref_GR)
	ol1_genes = unique(ol1@from)
	return(genes_GR[ol1_genes])
}

GENES_dm6_OUTDOM = myfindOverlaps(gene_dm6_gr, H3K27me3_OUTDOM)
GENES_dm6_OUTBORDER = myfindOverlaps(gene_dm6_gr, H3K27me3_OUTBORDER)
GENES_dm6_INDOM = myfindOverlaps(gene_dm6_gr, H3K27me3_INDOM)
GENES_dm6_INBORDER = myfindOverlaps(gene_dm6_gr, H3K27me3_INBORDER)


List_genes_DOM = list(GENES_dm6_OUTDOM = names(GENES_dm6_OUTDOM), GENES_dm6_OUTBORDER = names(GENES_dm6_OUTBORDER), GENES_dm6_INBORDER = names(GENES_dm6_INBORDER), GENES_dm6_INDOM = names(GENES_dm6_INDOM))

comRownames = Reduce(intersect,List_genes_DOM)
List_genes_DOM = lapply(List_genes_DOM, function(GN){GN = GN[GN %ni% comRownames]})

List_genes_DOM$GENES_dm6_OUTDOM = List_genes_DOM$GENES_dm6_OUTDOM[List_genes_DOM$GENES_dm6_OUTDOM %ni% c(List_genes_DOM$GENES_dm6_OUTBORDER, List_genes_DOM$GENES_dm6_INDOM, List_genes_DOM$GENES_dm6_INBORDER)]
List_genes_DOM$GENES_dm6_OUTBORDER = List_genes_DOM$GENES_dm6_OUTBORDER[List_genes_DOM$GENES_dm6_OUTBORDER %ni% c(List_genes_DOM$GENES_dm6_OUTDOM, List_genes_DOM$GENES_dm6_INDOM, List_genes_DOM$GENES_dm6_INBORDER)]
List_genes_DOM$GENES_dm6_INDOM = List_genes_DOM$GENES_dm6_INDOM[List_genes_DOM$GENES_dm6_INDOM %ni% c(List_genes_DOM$GENES_dm6_OUTDOM, List_genes_DOM$GENES_dm6_OUTBORDER, List_genes_DOM$GENES_dm6_INBORDER)]
List_genes_DOM$GENES_dm6_INBORDER = List_genes_DOM$GENES_dm6_INBORDER[List_genes_DOM$GENES_dm6_INBORDER %ni% c(List_genes_DOM$GENES_dm6_OUTDOM, List_genes_DOM$GENES_dm6_OUTBORDER, List_genes_DOM$GENES_dm6_INDOM)]


#  str(List_genes_DOM)
# List of 4
#  $ GENES_dm6_OUTDOM   : chr [1:5179] "FBgn0000017.1" "FBgn0000018.1" "FBgn0000043.1" "FBgn0000052.1" ...
#  $ GENES_dm6_OUTBORDER: chr [1:2893] "FBgn0000003.1" "FBgn0000032.1" "FBgn0000053.1" "FBgn0000064.1" ...
#  $ GENES_dm6_INBORDER : chr [1:3962] "FBgn0000024.1" "FBgn0000028.1" "FBgn0000036.1" "FBgn0000037.1" ...
#  $ GENES_dm6_INDOM    : chr [1:4948] "FBgn0000014.1" "FBgn0000015.1" "FBgn0000022.1" "FBgn0000038.1" ...


saveRDS(List_genes_DOM, paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/List_genes_DOM.RDS"))
List_genes_DOM = readRDS(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/List_genes_DOM.RDS"))

H3K27me3_CTRL_DOMAINS_borders_gapsup1500bp_sup1500bp_filledInputGaps
####################################################################################################################################################################
## END
####################################################################################################################################################################
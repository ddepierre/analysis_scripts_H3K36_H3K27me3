######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# SCRIPT PLOT HEATMAP FIGURE 1 D E HEATMAP

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
library("seqplots")
library("dplyr")
library("devtools")

## define workdir

workdir = ""
######################################################################################################################################################
### LOAD FUNCTION
######################################################################################################################################################


###### seqPlotSDoutliers
seqPlotSDoutliers <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,sd=10,err=T,type="pf",smooth=T,spar=0.35,gnme="dm6",ignore.strand=F, colvec = c("black","firebrick2")){
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(rtracklayer::export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }

  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
      gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
      gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
      means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
        sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
        qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      stderror[is.na(stderror)] <- 0
      conint[is.na(conint)] <- 0
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
    }

  }
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', ylim=ylim, ylab='Signal', main = paste0("Plot profile \n",sd," SD Removed"), keepratio = F,error.estimates = err, colvec = colvec, pointsize = 16, legend_ext = T)
  plotAverage(gpsa,xlab='Relative position [bp]', ylim=ylim, ylab='Signal', main = paste0("Plot profile \n",sd," SD Removed"), keepratio = F,error.estimates = err, colvec = colvec, pointsize = 16, legend=F)
}


## COLOR CODE 
# H3K27me3 #820002
# H3K36me2 #6ef4aa
# H3K36me3 #0d7c1c
#####################################################################################-
### LOAD BIGWIG
#####################################################################################-

tmp <- create(paste0(workdir,"FIGURES/FIGURE_1_B_BORDERDOM_AVG_RPOFILE/TMPgetPlotSetArray"))
tmp <- paste0(workdir,"FIGURES/FIGURE_1_B_BORDERDOM_AVG_RPOFILE/TMPgetPlotSetArray")

H3K36me2_RPGC = paste0(workdir, "DATA_PROCESSING/PROCESSED_BIGWIG_FILES/H3K36me2_chipseq_filt_sort_RPGC.bw")
H3K36me3_RPGC = paste0(workdir, "DATA_PROCESSING/PROCESSED_BIGWIG_FILES/H3K36me3_chipseq_filt_sort_RPGC.bw")
K27C_RPGC = paste0(workdir, "DATA_PROCESSING/PROCESSED_BIGWIG_FILES/K27C_trimmed_filt_sort_RPGC.bw")

################################################################################################################################################################################
# ref
gene_dm6_gr <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
seqlevelsStyle(gene_dm6_gr) <- "Ensembl"
names(gene_dm6_gr) = paste0(names(gene_dm6_gr),".1")
gene_dm6_gr <- gene_dm6_gr[seqnames(gene_dm6_gr) %in% c("2L", "2R", "3L", "3R", "4",  "X",  "Y")]
tss_gene_dm6_gr <- resize(gene_dm6_gr,1,"start")

H3K27me3_CTRL_DOMAINS_borders = rtracklayer::import.bed(paste0(workdir, "DATA_PROCESSING/NORMR_DETECTION/DOMAINS_H3K27me3/H3K27me3_CTRL_DOMAINS_borders_gapsup1500bp_sup1500bp_filledInputGaps.bed"))
GR_list_toPlot = c("H3K27me3_CTRL_DOMAINS_borders")

for(GR in GR_list_toPlot){
  pdf(paste0(workdir,"FIGURES/FIGURE_1_B_BORDERDOM_AVG_RPOFILE/AVERAGE_PROFILE_K27C_H3K36me2_H3K36me3_" ,GR,".pdf"))
  seqPlotSDoutliers(c(K27C_RPGC, H3K36me3_RPGC, H3K36me2_RPGC),tmp,GR,c(0,3),c(5000,5000),20L, colvec = c("#820002", "#0d7c1c", "#6ef4aa"))
  dev.off()
}


## END
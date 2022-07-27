######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# SCRIPT BOXPLOT SUPP FIGURE 6 B C 
# This Boxplot compares eigen values computed on HiC data, reflecting A/B compartment, in WT, MES4 KD and HypB KD for genes regulated by Mes-4 or HypB (i.e. genes with H3K27me3 spreading) 

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
library("ggridges")
library("ggthemes")
library("ggdark")
'%ni%' = Negate('%in%')


## define workdir
workdir = "/path/to/analysis_scripts_H3K36_H3K27me3/"



######################################################################################################################################################
### LOAD FUNCTION
######################################################################################################################################################


        # Plot Theme ofr boxplot
            my_theme = function (base_size = 11, base_family = 'Helvetica')
            {
            (theme_foundation(base_size = base_size, base_family = base_family) +
                theme(  line = element_line(colour = 'black'),
                        rect = element_rect(fill = 'white',linetype = 0, colour = NA), 
                        text = element_text(colour = 'black'),
                        axis.title = element_text(),
                        axis.text = element_text(),
                        axis.ticks = element_line(),
                        axis.line = element_line(size=1),
                        legend.background = element_rect(),
                        legend.position = 'bottom',
                        legend.direction = 'horizontal',
                        legend.box = 'horizontal',
                        legend.title = element_text(vjust=0.8),
                        panel.grid = element_line(colour = NULL),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        plot.title = element_text(hjust = 0, size=rel(1.5), face = 'bold'),
                        plot.tag = element_text(size=rel(1.8), face = 'bold'),
                        plot.subtitle = element_text(hjust = 0.2, size=rel(1.2), face = 'italic'),
                        plot.margin = unit(c(1,1, 1, 1), 'lines'),
                        plot.background = element_rect(linetype=1,colour='white',size=1.2),
                        strip.background = element_rect()))
            }
            


######################################################################################################################################################
### LOAD DATA
######################################################################################################################################################
### LOAD QUANTIF AND SCORES
QUANTIFandSCORES_LIST = readRDS(paste0(workdir, "DATA_PROCESSING/QUANTIFandSCORES_LIST.RDS"))
str(QUANTIFandSCORES_LIST)

### EIGEN VALUES
# Ramirez_compartment_res5000_S2_dm6 = readRDS(paste0(workdir, "DATA_PROCESSING/COMPARTMENT_AB/Ramirez_WT_AB_comp_res5000_S2_dm6_bed.RDS"))
# Ramirez_MES4_KD_AB_comp_res5000_S2_dm6 = readRDS(paste0(workdir, "DATA_PROCESSING/COMPARTMENT_AB/Ramirez_MES4_KD_AB_comp_res5000_S2_dm6_bed.RDS"))
eigenvector_HYPB_KD_mergeRep_ABcomp_res5000 = readRDS(paste0(workdir, "DATA_PROCESSING/COMPARTMENT_AB/eigenvector_HYPB_KD_mergeRep_inter_30_compartment_AB_comp_res5000_S2_dm6_bed.RDS"))
eigenvector_luc_KD_rep1_ABcomp_res5000 = readRDS(paste0(workdir, "DATA_PROCESSING/COMPARTMENT_AB/eigenvector_luc_KD_rep1_inter_30_compartment_AB_comp_res5000_S2_dm6_bed.RDS"))
eigenvector_MES4_KD_mergeRep_ABcomp_res5000 = readRDS(paste0(workdir, "DATA_PROCESSING/COMPARTMENT_AB/eigenvector_MES4_KD_mergeRep_inter_30_compartment_AB_comp_res5000_S2_dm6_bed.RDS"))



# Load ref genes Granges file
gene_dm6_gr = readRDS(paste0(workdir, "DATA_PROCESSING/REF_GENOME/TxDb.Dmelanogaster.Ensembl.dm6.RDS"))
names(gene_dm6_gr) = paste0(names(gene_dm6_gr), ".1")
seqlevelsStyle(gene_dm6_gr) = "UCSC"

#####################################################################################################################################################
### PREPARE DATA TO PLOT
######################################################################################################################################################
## OVERLAP GENES WITH COMPARTMENTS
## WT
gene_dm6_gr_TSS = resize(gene_dm6_gr,1, fix="start")
gene_dm6_gr_TSS_comp_ovlp = findOverlaps(gene_dm6_gr_TSS, eigenvector_luc_KD_rep1_ABcomp_res5000)
COMPARTMENT_SCORE_gene_dm6_gr_TSS = eigenvector_luc_KD_rep1_ABcomp_res5000[gene_dm6_gr_TSS_comp_ovlp@to]$score
names(COMPARTMENT_SCORE_gene_dm6_gr_TSS) = paste0(gene_dm6_gr_TSS[gene_dm6_gr_TSS_comp_ovlp@from]$gene_id, ".1")

## Mes-4 KD
gene_dm6_gr_TSS = resize(gene_dm6_gr,1, fix="start")
gene_dm6_gr_TSS_comp_MES4_KD_ovlp = findOverlaps(gene_dm6_gr_TSS, eigenvector_MES4_KD_mergeRep_ABcomp_res5000)
COMPARTMENT_SCORE_MES4_KD_gene_dm6_gr_TSS = eigenvector_MES4_KD_mergeRep_ABcomp_res5000[gene_dm6_gr_TSS_comp_MES4_KD_ovlp@to]$score
names(COMPARTMENT_SCORE_MES4_KD_gene_dm6_gr_TSS) = paste0(gene_dm6_gr_TSS[gene_dm6_gr_TSS_comp_MES4_KD_ovlp@from]$gene_id, ".1")

## HypB KD
gene_dm6_gr_TSS = resize(gene_dm6_gr,1, fix="start")
gene_dm6_gr_TSS_comp_HYPB_KD_ovlp = findOverlaps(gene_dm6_gr_TSS, eigenvector_HYPB_KD_mergeRep_ABcomp_res5000)
COMPARTMENT_SCORE_HYPB_KD_gene_dm6_gr_TSS = eigenvector_HYPB_KD_mergeRep_ABcomp_res5000[gene_dm6_gr_TSS_comp_HYPB_KD_ovlp@to]$score
names(COMPARTMENT_SCORE_HYPB_KD_gene_dm6_gr_TSS) = paste0(gene_dm6_gr_TSS[gene_dm6_gr_TSS_comp_HYPB_KD_ovlp@from]$gene_id, ".1")


# GET ZSCORE_K27M_K27H
ZSCORE_K27M_K27H = QUANTIFandSCORES_LIST$ZSCORE_K27M_K27H
Nsplit = 5
ZSCORE_K27M_K27H_GN_splitted = split(names(ZSCORE_K27M_K27H), ceiling(seq_along(names(ZSCORE_K27M_K27H))/ceiling(length(names(ZSCORE_K27M_K27H))/Nsplit)))
ZSCORE_K27M_K27H_VALUE_splitted = split(ZSCORE_K27M_K27H, ceiling(seq_along(ZSCORE_K27M_K27H)/ceiling(length(ZSCORE_K27M_K27H)/Nsplit)))

####################################################################################################################################################################################################################
### BOXPLOT WT vs Mes-4 KD eigen values
####################################################################################################################################################################################################################



filt_names_UP_K27 = list(
UP_K27_MES4_SPECIFIC = ZSCORE_K27M_K27H_GN_splitted[[1]],
UP_K27_HYPB_SPECIFIC = ZSCORE_K27M_K27H_GN_splitted[[5]],
Randoms	= sample(names(gene_dm6_gr), 1000)
)




        couplesTOplot = c(
                "UP_K27_MES4_SPECIFIC",
                "UP_K27_HYPB_SPECIFIC",
                "Randoms"
                )
 
        
        color.vec=c("#00b0f0", "#002060", "#a6a6a6")
        
        X_vals = "Class" # Metric Column name to use
        # NORM_VAL = "AllTSS_AllTSS"
        NORM_VAL = "Randoms"
        VALUEtoplot = "COMPARTMENT_SCORE"
        NAMEtoplot = paste0(VALUEtoplot, "_MES4", "_", paste(couplesTOplot, collapse="__"))
        YLIM = c(-25, 25)

                            plot.name = paste0("BOXPLOT_",NAMEtoplot)

							metricDataFrame_WT = data.frame(
								Individus = as.character(names(unlist(filt_names_UP_K27))),
								COMPARTMENT_SCORE = c(COMPARTMENT_SCORE_gene_dm6_gr_TSS[filt_names_UP_K27$UP_K27_MES4_SPECIFIC],
													COMPARTMENT_SCORE_gene_dm6_gr_TSS[filt_names_UP_K27$UP_K27_HYPB_SPECIFIC],
													COMPARTMENT_SCORE_gene_dm6_gr_TSS[filt_names_UP_K27$Randoms])

							)
							metricDataFrame_WT$Individus = as.character(metricDataFrame_WT$Individus)

							metricDataFrame_KD = data.frame(
								Individus = as.character(names(unlist(filt_names_UP_K27))),
								COMPARTMENT_SCORE = c(COMPARTMENT_SCORE_MES4_KD_gene_dm6_gr_TSS[filt_names_UP_K27$UP_K27_MES4_SPECIFIC],
													COMPARTMENT_SCORE_MES4_KD_gene_dm6_gr_TSS[filt_names_UP_K27$UP_K27_HYPB_SPECIFIC],
													COMPARTMENT_SCORE_MES4_KD_gene_dm6_gr_TSS[filt_names_UP_K27$Randoms])

							)
							metricDataFrame_KD$Individus = as.character(metricDataFrame_KD$Individus)

                                      # FILTER metricDataFrame_WT
                            dim(metricDataFrame_WT)
                            metricDataFrame_WT$Class = gsub('[0-9]{1,4}$', '', metricDataFrame_WT$Individus)
                            metricDataFrame_WT$Class = as.factor(metricDataFrame_WT$Class)
                            metricDataFrame_WT = metricDataFrame_WT[metricDataFrame_WT$Class %in% couplesTOplot,] 
                            metricDataFrame_WT$Class = factor(metricDataFrame_WT$Class, levels=couplesTOplot)
                            dim(metricDataFrame_WT)

                                      # FILTER metricDataFrame_KD
                            dim(metricDataFrame_KD)
                            metricDataFrame_KD$Class = gsub('[0-9]{1,4}$', '', metricDataFrame_KD$Individus)
                            metricDataFrame_KD$Class = as.factor(metricDataFrame_KD$Class)
                            metricDataFrame_KD = metricDataFrame_KD[metricDataFrame_KD$Class %in% couplesTOplot,] 
                            metricDataFrame_KD$Class = factor(metricDataFrame_KD$Class, levels=couplesTOplot)
                            dim(metricDataFrame_KD)

                                # FILTER metricDataFrame_KD

                            metricDataFrame = rbind(metricDataFrame_WT, metricDataFrame_KD)
                            dim(metricDataFrame)

                            metricDataFrame$COND = c(rep("WT",nrow(metricDataFrame)/2), rep("MES4_KD", nrow(metricDataFrame)/2))
                            metricDataFrame$COND = as.factor(metricDataFrame$COND)
                            metricDataFrame$COND = factor(metricDataFrame$COND, levels = c("WT", "MES4_KD"))

                            # define wilcox complist
                            comp_list = list(c(levels(metricDataFrame$Class)[1], levels(metricDataFrame$Class)[2]), 
                                             c(levels(metricDataFrame$Class)[1], levels(metricDataFrame$Class)[3]),
                                             c(levels(metricDataFrame$Class)[2], levels(metricDataFrame$Class)[3]),
                                             c("MES4_KD", "WT"))

                            ## REDEF SCORES AS NUMERIC      
                            metricDataFrame$Individus=NULL
                                    
                            ### NORMALIZED VALUES
                            avgCTRL_WT = mean(metricDataFrame[metricDataFrame$Class %in% NORM_VAL & metricDataFrame$COND %in% "WT",VALUEtoplot], na.rm=T)
                            avgCTRL_KD = mean(metricDataFrame[metricDataFrame$Class %in% NORM_VAL & metricDataFrame$COND %in% "MES4_KD",VALUEtoplot], na.rm=T)
                            normFact = avgCTRL_WT/avgCTRL_KD
                            
                            metricDataFrame_norm = metricDataFrame
                            metricDataFrame_norm[metricDataFrame$COND %in% "MES4_KD",VALUEtoplot] = metricDataFrame_norm[metricDataFrame$COND %in% "MES4_KD",VALUEtoplot]*normFact

                            ##############
                            ## PLOT BOXPLOT
                            ##############
                            
                            boxplot_ttest_norm = ggplot(metricDataFrame_norm, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) + 
                            scale_fill_manual(values=color.vec) + scale_color_manual(values=c("#525252", "#b52828")) +
                            labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("T.test + Norm on : ", NORM_VAL)) + stat_compare_means(aes(group = COND), method = "t.test", label = "p.format") + my_theme()
                                
                            boxplotYLIM_norm = ggplot(metricDataFrame_norm, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) + 
                            scale_fill_manual(values=color.vec) + ylim(YLIM) + scale_color_manual(values=c("#525252", "#b52828")) +
                            labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("FORCED YLIM + Norm on : ", NORM_VAL)) + my_theme()


                            print(paste0(workdir, "FIGURES/FIGURE_5_B_&_Supp_6_BC_BOXPLOT_EIGEN_VALUE/",plot.name,".pdf"))
                            pdf(file=paste0(workdir, "FIGURES/FIGURE_5_B_&_Supp_6_BC_BOXPLOT_EIGEN_VALUE/",plot.name,".pdf"))
                            print(boxplot_ttest_norm)
                            print(boxplotYLIM_norm)
                            dev.off()


####################################################################################################################################################################################################################
### BOXPLOT WT vs HypB KD eigen values
####################################################################################################################################################################################################################



filt_names_UP_K27 = list(
UP_K27_MES4_SPECIFIC = ZSCORE_K27M_K27H_GN_splitted[[1]],
UP_K27_HYPB_SPECIFIC = ZSCORE_K27M_K27H_GN_splitted[[5]],
Randoms	= sample(names(gene_dm6_gr), 1000)
)




        couplesTOplot = c(
                "UP_K27_MES4_SPECIFIC",
                "UP_K27_HYPB_SPECIFIC",
                "Randoms"
                )
 
        
        color.vec=c("#00b0f0", "#002060", "#a6a6a6")
        
        X_vals = "Class" # Metric Column name to use
        # NORM_VAL = "AllTSS_AllTSS"
        NORM_VAL = "Randoms"
        VALUEtoplot = "COMPARTMENT_SCORE"
        NAMEtoplot = paste0(VALUEtoplot, "_HYPB", "_", paste(couplesTOplot, collapse="__"))
        YLIM = c(-25, 25)

                            plot.name = paste0("BOXPLOT_",NAMEtoplot)

							metricDataFrame_WT = data.frame(
								Individus = as.character(names(unlist(filt_names_UP_K27))),
								COMPARTMENT_SCORE = c(COMPARTMENT_SCORE_gene_dm6_gr_TSS[filt_names_UP_K27$UP_K27_MES4_SPECIFIC],
													COMPARTMENT_SCORE_gene_dm6_gr_TSS[filt_names_UP_K27$UP_K27_HYPB_SPECIFIC],
													COMPARTMENT_SCORE_gene_dm6_gr_TSS[filt_names_UP_K27$Randoms])

							)
							metricDataFrame_WT$Individus = as.character(metricDataFrame_WT$Individus)

							metricDataFrame_KD = data.frame(
								Individus = as.character(names(unlist(filt_names_UP_K27))),
								COMPARTMENT_SCORE = c(COMPARTMENT_SCORE_HYPB_KD_gene_dm6_gr_TSS[filt_names_UP_K27$UP_K27_MES4_SPECIFIC],
													COMPARTMENT_SCORE_HYPB_KD_gene_dm6_gr_TSS[filt_names_UP_K27$UP_K27_HYPB_SPECIFIC],
													COMPARTMENT_SCORE_HYPB_KD_gene_dm6_gr_TSS[filt_names_UP_K27$Randoms])

							)
							metricDataFrame_KD$Individus = as.character(metricDataFrame_KD$Individus)

                                      # FILTER metricDataFrame_WT
                            dim(metricDataFrame_WT)
                            metricDataFrame_WT$Class = gsub('[0-9]{1,4}$', '', metricDataFrame_WT$Individus)
                            metricDataFrame_WT$Class = as.factor(metricDataFrame_WT$Class)
                            metricDataFrame_WT = metricDataFrame_WT[metricDataFrame_WT$Class %in% couplesTOplot,] 
                            metricDataFrame_WT$Class = factor(metricDataFrame_WT$Class, levels=couplesTOplot)
                            dim(metricDataFrame_WT)

                                      # FILTER metricDataFrame_KD
                            dim(metricDataFrame_KD)
                            metricDataFrame_KD$Class = gsub('[0-9]{1,4}$', '', metricDataFrame_KD$Individus)
                            metricDataFrame_KD$Class = as.factor(metricDataFrame_KD$Class)
                            metricDataFrame_KD = metricDataFrame_KD[metricDataFrame_KD$Class %in% couplesTOplot,] 
                            metricDataFrame_KD$Class = factor(metricDataFrame_KD$Class, levels=couplesTOplot)
                            dim(metricDataFrame_KD)

                                # FILTER metricDataFrame_KD

                            metricDataFrame = rbind(metricDataFrame_WT, metricDataFrame_KD)
                            dim(metricDataFrame)

                            metricDataFrame$COND = c(rep("WT",nrow(metricDataFrame)/2), rep("HYPB_KD", nrow(metricDataFrame)/2))
                            metricDataFrame$COND = as.factor(metricDataFrame$COND)
                            metricDataFrame$COND = factor(metricDataFrame$COND, levels = c("WT", "HYPB_KD"))

                            # define wilcox complist
                            comp_list = list(c(levels(metricDataFrame$Class)[1], levels(metricDataFrame$Class)[2]), 
                                             c(levels(metricDataFrame$Class)[1], levels(metricDataFrame$Class)[3]),
                                             c(levels(metricDataFrame$Class)[2], levels(metricDataFrame$Class)[3]),
                                             c("HYPB_KD", "WT"))

                            ## REDEF SCORES AS NUMERIC      
                            metricDataFrame$Individus=NULL
                                    
                            ### NORMALIZED VALUES
                            avgCTRL_WT = mean(metricDataFrame[metricDataFrame$Class %in% NORM_VAL & metricDataFrame$COND %in% "WT",VALUEtoplot], na.rm=T)
                            avgCTRL_KD = mean(metricDataFrame[metricDataFrame$Class %in% NORM_VAL & metricDataFrame$COND %in% "HYPB_KD",VALUEtoplot], na.rm=T)
                            normFact = avgCTRL_WT/avgCTRL_KD
                            
                            metricDataFrame_norm = metricDataFrame
                            metricDataFrame_norm[metricDataFrame$COND %in% "HYPB_KD",VALUEtoplot] = metricDataFrame_norm[metricDataFrame$COND %in% "HYPB_KD",VALUEtoplot]*normFact

                            ##############
                            ## PLOT BOXPLOT
                            ##############
                            
                            # boxplot = ggplot(metricDataFrame, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) +  
                            # scale_fill_manual(values=color.vec) + scale_color_manual(values=c("#525252", "#b52828")) +
                            # labs(title=plot.name, x=X_vals, y=print(VALUEtoplot)) + my_theme()

                            # boxplot_wilcox = ggplot(metricDataFrame, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) + 
                            # scale_fill_manual(values=color.vec) + scale_color_manual(values=c("#525252", "#b52828")) +
                            # labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("WILCOXON")) + stat_compare_means(aes(group = COND), method = "wilcox.test", label = "p.format") + my_theme()

                            # boxplot_ttest = ggplot(metricDataFrame, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) + 
                            # scale_fill_manual(values=color.vec) + scale_color_manual(values=c("#525252", "#b52828")) +
                            # labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("T.test")) + stat_compare_means(aes(group = COND), method = "t.test", label = "p.format") + my_theme()

                            # boxplotYLIM = ggplot(metricDataFrame, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) + 
                            # scale_fill_manual(values=color.vec) + ylim(YLIM) + scale_color_manual(values=c("#525252", "#b52828")) +
                            # labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("FORCED YLIM")) + my_theme()

                            # # metricDataFrame_norm 
                            # boxplot_norm = ggplot(metricDataFrame_norm, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) +  
                            # scale_fill_manual(values=color.vec) + scale_color_manual(values=c("#525252", "#b52828")) +
                            # labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("Norm on : ", NORM_VAL)) + my_theme()

                            # boxplot_wilcox_norm = ggplot(metricDataFrame_norm, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) + 
                            # scale_fill_manual(values=color.vec) + scale_color_manual(values=c("#525252", "#b52828")) +
                            # labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("WILCOXON + Norm on : ", NORM_VAL)) + stat_compare_means(aes(group = COND), method = "wilcox.test", label = "p.format") + my_theme()

                            boxplot_ttest_norm = ggplot(metricDataFrame_norm, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) + 
                            scale_fill_manual(values=color.vec) + scale_color_manual(values=c("#525252", "#b52828")) +
                            labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("T.test + Norm on : ", NORM_VAL)) + stat_compare_means(aes(group = COND), method = "t.test", label = "p.format") + my_theme()
                                
                            boxplotYLIM_norm = ggplot(metricDataFrame_norm, aes(  x = get(X_vals), y=  get(VALUEtoplot), fill=get(X_vals), color=COND))  + geom_boxplot(lwd=2) + 
                            scale_fill_manual(values=color.vec) + ylim(YLIM) + scale_color_manual(values=c("#525252", "#b52828")) +
                            labs(title=plot.name, x=X_vals, y=print(VALUEtoplot), subtitle = paste0("FORCED YLIM + Norm on : ", NORM_VAL)) + my_theme()


                            print(paste0(workdir, "FIGURES/FIGURE_5_B_&_Supp_6_BC_BOXPLOT_EIGEN_VALUE/",plot.name,".pdf"))
                            pdf(file=paste0(workdir, "FIGURES/FIGURE_5_B_&_Supp_6_BC_BOXPLOT_EIGEN_VALUE/",plot.name,".pdf"))
                            # print(boxplot)
                            # print(boxplot_wilcox)
                            # print(boxplot_ttest)
                            # print(boxplotYLIM)
                            # print(boxplot_norm)
                            # print(boxplot_wilcox_norm)
                            print(boxplot_ttest_norm)
                            print(boxplotYLIM_norm)
                            dev.off()






#end

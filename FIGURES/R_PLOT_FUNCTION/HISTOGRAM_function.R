######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

## FUNCTION TO COMPARE DISTRIBUTION of a same feature between 2 list of IDs
## input : vec = named numeric vector, vec_GN = subset of vec names

print('USAGE : plotHisto_listGN(vec = vec, nameVec = "nameVec", vec_GN = vec_GN, nameVec_GN = "nameVec_GN", xlim = c(0, 30000), ylim = c(0, 0.0002), fillingcol = c("#3d3d3d", "#820002"),
outDir = paste0(workdir, "path/to/out/"), info = NULL)')


######################################################################################################################################################
######################################################################################################################################################
# R version 3.4.2 (2017-09-28) -- "Short Summer"
# Copyright (C) 2017 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)



##################################################################################
## LIBRARY
##################################################################################

library("data.table")
library("ggplot2")


##################################################################################
## FUNCTION
##################################################################################

plotHisto_listGN = function(vec = vec, nameVec = "nameVec", vec_GN = vec_GN, nameVec_GN = "nameVec_GN", xlim = c(0, 30000), ylim = c(0, 0.0002), binW = 0.5, alphaTransp = 0.6, fillingcol = c("#3d3d3d", "#820002"),
outDir = paste0(workdir, "path/to/out/"), info = NULL){
  data2GG = data.table(vec)
  rownames(data2GG) = names(vec)
  data2GG$feat1 = "A"
  data2GG$feat1[rownames(data2GG) %in% vec_GN] = "B"
  data2GG$feat1 = as.factor(data2GG$feat1)
  pdf(paste0(outDir, "HISTOGRAM_",nameVec ,"_Spllited_", nameVec_GN,"_",info,".pdf"), width = 16, height = 10 )
  p_hist = ggplot(data=data2GG, aes(data2GG$vec, fill = data2GG$feat1, stat(density))) + xlim(xlim) + ylim(ylim) + xlab(nameVec) +
  geom_histogram(binwidth=binW,alpha=alphaTransp, position="identity")  +
  scale_fill_manual(values = fillingcol,
    name=nameVec_GN,
    labels=c("NO", "YES")) +
    theme_classic(base_line_size = 1, base_rect_size=1) +
    geom_vline(data=data2GG, aes(xintercept=0), linetype="dashed", size=1, colour="#3d3d3d")
  p_hist2 = ggplot(data=data2GG, aes(data2GG$vec, fill = data2GG$feat1, stat(density))) + xlim(xlim) + ylim(ylim) + xlab(nameVec) +
  geom_histogram(binwidth=binW,alpha=1, position="dodge")  +
  scale_fill_manual(values = fillingcol,
    name=nameVec_GN,
    labels=c("NO", "YES")) +
    theme_classic(base_line_size = 1, base_rect_size=1)
    print(p_hist)
    print(p_hist2)

    dev.off()
}



#end

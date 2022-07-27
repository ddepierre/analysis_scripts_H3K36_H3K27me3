######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

# From 2 lists of IDs list1 and list2
# 1/ get intersection matrices of each 1st list elements in 2nd list elements => FUNCTION "fisher_namesList"
# out put matrices are intersection effectives & intersection percentages & fisher exact test intersection testing (odds ratio and pvalues) 
# 2/ Get a graphic visualisation of intersection by plotting a heatmap from Fisher odss ratio or pvalues => FUNCTION "Hmap_pval"

######################################################################################################################################################
######################################################################################################################################################
# R version 3.4.2 (2017-09-28) -- "Short Summer"
# Copyright (C) 2017 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

######################################################################################################################################################
######################################################################################################################################################

library("pheatmap")



fisher_namesList=function(list1, list2, total, test.side = "greater")
 {
 require(GenomicRanges)
 ## create matrices to fill
 m_pval=matrix(nrow=length(list1),ncol=length(list2))
 m_or=matrix(nrow=length(list1),ncol=length(list2))
 m_eff=matrix(nrow=length(list1),ncol=length(list2))
 m_pct_list1=matrix(nrow=length(list1),ncol=length(list2))
 m_pct_list2=matrix(nrow=length(list1),ncol=length(list2))

 colnames(m_pval)=names(list2)
 rownames(m_pval)=names(list1)
 colnames(m_or)=names(list2)
 rownames(m_or)=names(list1)
 colnames(m_eff)=names(list2)
 rownames(m_eff)=names(list1)
 colnames(m_pct_list1)=names(list2)
 rownames(m_pct_list1)=names(list1)
 colnames(m_pct_list2)=names(list2)
 rownames(m_pct_list2)=names(list1)

 ## fill matrices with either fisher pval, fisher odds, effectives or percentage
 for (i in 1:length(list1)){
   liste1=list1[[i]]
    for (j in 1:length(list2)){
      liste2=list2[[j]]
      inter=length(Reduce(intersect, list(liste1,liste2)))
      m_pval[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative=test.side)$p.value,2)
      m_or[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative=test.side)$estimate,2)
      m_eff[i,j]=inter
      m_pct_list1[i,j] = inter/length(liste1)
      m_pct_list2[i,j] = inter/length(liste2)

    }
 }

## get list of matrices to return
 l=list()
 l[[1]]=m_pval
 l[[2]]=m_or
 l[[3]]=m_eff
 l[[4]]=m_pct_list1
 l[[5]]=m_pct_list2
 names(l)=c("p.value", "odds.ratio", "effectives", "pct_list1", "pct_list2")
 return(l)
 }


########################################################################################################################################################################
###    FROM MATRICES RETURN FROM fisher_namesList(), PLOT HEATMAP
########################################################################################################################################################################
 Hmap_pval=function(mat, values_ordonnée, values_abscisse, graph,bool, threshold=1e-8, odds_scale = c(0,2)){
   if (bool==T){pdf(graph)} # should we plot in pdf -> yes if bool = T
   # Prepare pval matrix
   mat_pv=mat$p.value
   mat_pv=mat_pv+10^-300 # to avoid pval = 0
   # ADDED on 23/10/2019
   mat_pv[mat_pv > 0.05] = 1 # set pval to 1 if >0.05 (basis sup threshold)
   if (threshold>1 | threshold<=0) (stop("threshold must be in ]0,1]"))
   mat_pv[mat_pv < threshold] = threshold # set pval to threshold if threshold (basis inferior threshold)
   mat_pv=log10(mat_pv) # transform pval to log10


##########################
####   PLOT PVALUE   #####
##########################
   par(mfrow=c(2,1))
   logt=log10(threshold)
   pval_info = matrix(c(1, threshold))
   rownames(pval_info) = c("white", "blue")
   textplot(pval_info)
   textplot(mat$p.value)
   title("PVALUES FISHER EXACT TEST")

   par(mfrow=c(1,1))

  v_breaks = c(seq(min(mat_pv), 0.1, length.out=25), seq(0.09, max(mat_pv), length.out=25))
  v_cols = c(colorRampPalette(c("blue", "white"))(24), "white", colorRampPalette(c("white", "white"))(24))
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = mat$p.value, fontsize_number = 18, main="PVALUES FISHER EXACT TEST")
  v_breaks = c(seq(log10(threshold), 0.1, length.out=25), seq(0.09, max(mat_pv), length.out=25))
  v_cols = c(colorRampPalette(c("blue", "white"))(24), "white", colorRampPalette(c("white", "white"))(24))
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, main="PVALUES FISHER EXACT TEST")
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = mat$p.value, fontsize_number = 18, main="PVALUES FISHER EXACT TEST")
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = mat$effectives, fontsize_number = 18, main="PVALUES FISHER EXACT TEST - EFFECTIVES")
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = round(mat$pct_list1,2), fontsize_number = 18, main="% of list2 in list1 / i.e. column feat. in row feat.")
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = round(mat$pct_list2,2), fontsize_number = 18, main="% of list1 in list2 / i.e. row feat. in comumn feat.")


#############################
####   PLOT ODDS RATIO   ####
#############################
    textplot(mat$odds.ratio) # PLOT ODDS RATIO TABLE
    title("ODDS RATIO FISHER EXACT TEST")

    plot.new()
    ## RAW ODDS RATIO
    M = mat$odds.ratio
    v_breaks = c(seq(min(M), 0.9, length.out=25), seq(1.1, max(M), length.out=25))
    v_cols = c(colorRampPalette(c("yellow", "white"))(24), "white", colorRampPalette(c("white", "darkblue"))(24))
    pheatmap(M, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, main="ODDS RATIO FISHER EXACT TEST")
    pheatmap(M, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, display_numbers = mat$odds.ratio, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST")
    pheatmap(M, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, display_numbers = mat$p.value, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST - pval numbers displayed")
    pheatmap(M, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = mat$effectives, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST - EFFECTIVES")


    ## LIMITED ODDS RATIO to borders defined with "odds_scale" param
    M_v2 = M
    M_v2[M_v2 > odds_scale[2]] = odds_scale[2]
    v_breaks = c(seq(odds_scale[1], 0.8, length.out=25), seq(1.2, odds_scale[2], length.out=25))
    v_cols = c(colorRampPalette(c("yellow", "white"))(24), "white", colorRampPalette(c("white", "darkblue"))(24))
    pheatmap(M_v2, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, main="ODDS RATIO FISHER EXACT TEST")
    pheatmap(M_v2, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, display_numbers = mat$odds.ratio, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST")
    pheatmap(M_v2, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, display_numbers = mat$p.value, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST - pval numbers displayed")


    if(is.character(values_ordonnée[[1]])){ #if the split of genes nmaes is not done on quantitatif value, plot effectif of each feature
      barplot(unlist(lapply(values_ordonnée, length)), cex.names=0.5)
    }else{
      boxplot(values_ordonnée, main=deparse(substitute(values_ordonnée)), outline=F)
      barplot(unlist(values_ordonnée), main=deparse(substitute(values_ordonnée)))
      plot(density(unlist(values_ordonnée)), main=deparse(substitute(values_ordonnée)))
      hist(unlist(values_ordonnée), 50, main=deparse(substitute(values_ordonnée)))
    }
    if(is.character(values_abscisse[[1]])){ #if the split of genes nmaes is not done on quantitatif value, plot effectif of each feature
      barplot(unlist(lapply(values_abscisse, length)), cex.names=0.5)
    }else{
      boxplot(values_abscisse, main=deparse(substitute(values_abscisse)), outline=F)
      barplot(unlist(values_abscisse), main=deparse(substitute(values_abscisse)))
      plot(density(unlist(values_abscisse)), main=deparse(substitute(values_abscisse)))
      hist(unlist(values_abscisse), 50, main=deparse(substitute(values_abscisse)))
    }
   if (bool==T){
     dev.off()
   }
 }


# end

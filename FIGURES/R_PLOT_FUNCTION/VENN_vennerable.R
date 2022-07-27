######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################


# FUNCTION TO PLOT PROPORTIONAL VENN DIRAGRAM FROM 2 or 3 names vector stored in a list

######################################################################################################################################################
######################################################################################################################################################

print('From list of 2 or 3 character vectors, plot venn + enrichment test fisher + SuperExactTest
				USAGE : plotVENNerable(data_to_venn = list(vec1, vec2, vec3), Nb_REF = length(univers), outdir = paste0(workdir, "path/to/plot/"))')



######################################################################################################################################################
## LIBRARY
######################################################################################################################################################
library(gplots)
library(Vennerable)
library(SuperExactTest)

######################################################################################################################################################
## FUNCTION
######################################################################################################################################################

plotVENNerable = function(data_to_venn, Nb_REF, outdir){
	vennNAME = paste(names(data_to_venn), collapse="_")
	Venn_ovlp = Venn(data_to_venn)
	color_Venn = VennThemes(compute.Venn(Venn_ovlp))
  #set1
  color_Venn[["Set"]][["Set1"]]$col = "#016000"
  color_Venn[["SetText"]][["Set1"]]$col = "#016000"
  color_Venn[["Face"]][["100"]]$fill = "#63bf61"
  color_Venn[["Face"]][["100-1"]]$fill = "#63bf61"
  #set2
  # color_Venn[["Set"]][["Set2"]]$col = "#6d0000"
  # color_Venn[["SetText"]][["Set2"]]$col = "#6d0000"
  # color_Venn[["Face"]][["010"]]$fill = "#bc5656"
  # color_Venn[["Face"]][["010-1"]]$fill = "#bc5656"
  #set2
  color_Venn[["Set"]][["Set2"]]$col = "#ad6a05"
  color_Venn[["SetText"]][["Set2"]]$col = "#ad6a05"
  color_Venn[["Face"]][["010"]]$fill = "#fc9d44"
  color_Venn[["Face"]][["010-1"]]$fill = "#fc9d44"
  #set3
  color_Venn[["Set"]][["Set3"]]$col = #232323"
  color_Venn[["SetText"]][["Set3"]]$col = #232323"
  color_Venn[["Face"]][["001"]]$fill = "#a3a3a3"
  color_Venn[["Face"]][["001-1"]]$fill = "#a3a3a3"
  #set1 n set2
  color_Venn[["Face"]][["110"]]$fill = "#b5a029"
  color_Venn[["Face"]][["110-1"]]$fill = "#b5a029"
  #set1 n set3
  color_Venn[["Face"]][["101"]]$fill = "#276825"
  color_Venn[["Face"]][["101-1"]]$fill = "#276825"
  #set2 n set3
	color_Venn[["Face"]][["011"]]$fill = "#84603f"
  color_Venn[["Face"]][["011-1"]]$fill = "#84603f"
  # #set2 n set3
	# color_Venn[["Face"]][["011"]]$fill = "#721e1e"
  # color_Venn[["Face"]][["011-1"]]$fill = "#721e1e"

	# saveRDS(Venn_ovlp, paste0(outdir, "OBJvenn_",vennNAME,".RDS"))
	#fisher exact test
	fishertest = fisher_namesList(data_to_venn,data_to_venn,Nb_REF)
	# Supertest for triple intersection
	ResSuperTest = supertest(data_to_venn, n =Nb_REF)
	SupetTestToPrint = as.matrix(ResSuperTest$P.value)
	colnames(SupetTestToPrint) = "intersect_pval"
	########### PLOT #########################
	pdf(paste0(outdir, "PLOTvenn_",vennNAME,".pdf"))
	plot(Venn_ovlp, gp = color_Venn, doWeights = TRUE)
	par(mfrow=c(2,1))
	textplot(fishertest$p.value,  valign="top")
	textplot(fishertest$odds.ratio)
	textplot(SupetTestToPrint)
	plot(Venn_ovlp, gp = color_Venn, doWeights = TRUE, show = list(SetLabels = FALSE, FaceText = ""))
	dev.off()
}

### Compute fisher exact est table odds and pval from list venn
fisher_namesList=function(list1, list2, total, effMax = F)
 {
 require(GenomicRanges)
 if(effMax == T){
   effmax_l1 = min(unlist(lapply(list1, length)))
   list1 = lapply(list1, function(vec){
   	sample(vec, effmax_l1)
   })
   effmax_l2 = min(unlist(lapply(list2, length)))
   list2 = lapply(list2, function(vec){
   	sample(vec, effmax_l2)
   })
 }
 m_pval=matrix(nrow=length(list1),ncol=length(list2))
 m_or=matrix(nrow=length(list1),ncol=length(list2))
 for (i in 1:length(list1)){
   liste1=list1[[i]]
    for (j in 1:length(list2)){
      liste2=list2[[j]]
      inter=length(Reduce(intersect, list(liste1,liste2)))
      m_pval[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative="greater")$p.value,2)
      m_or[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative="greater")$estimate,2)
    }
 }
 colnames(m_pval)=names(list2)
 rownames(m_pval)=names(list1)
 colnames(m_or)=names(list2)
 rownames(m_or)=names(list1)
 l=list()
 l[[1]]=m_pval
 l[[2]]=m_or
 names(l)=c("p.value", "odds.ratio")
 return(l)
 }

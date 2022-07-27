######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################


# FUNCTION TO DO PCA AND DENDROGRAM CLUSTERING on PCs projection as distances
# Thanks to Raoul Raffel

######################################################################################################################################################
## LIBRARY
######################################################################################################################################################

invisible(require(dendsort))
invisible(require(FactoMineR))
invisible(require(dendextend))


######################################################################################################################################################
## FUNCTION
######################################################################################################################################################

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

do_acp_and_clustering <- function(
                            dat,
                            list.clust=list(c(1,2),c(1,3),c(2,3)),
                            list.cpa=list(c(1,2),c(1,3),c(2,3)),
                            cpa.title="ACP", dist.method="pearson",
                            choix="ind",
                            col.ind = "1"){

  invisible(require(FactoMineR))
  invisible(require(dendextend))


  old.par <- par()
  if(class(dat) %in% "list"){
    comRownames = Reduce(intersect,lapply(dat,names))
    dat = lapply(dat, function(feat1){feat1 = feat1[names(feat1) %in% comRownames]})
    data_df = matrix(NA,ncol = length(dat), nrow = length(dat[[1]]))
    rownames(data_df) = names(dat[[1]])
    colnames(data_df) = names(dat)
    for(i in 1:length(dat)){
      data_df[,i] = dat[[i]][rownames(data_df)]
    }
    dat = data_df
  }

  stopifnot(choix %in% c("var", "ind"))

  ncp <- max(c(unlist(list.clust), unlist(list.cpa)))
  ncp <- min(ncp, dim(dat)[2])

  dat.PCA=PCA(X=dat, graph=F, ncp=ncp)
  eig_title=paste(c("eigenvalue for",cpa.title),collapse=" ")
  eig <- dat.PCA$eig[1:ncp,"percentage of variance"]
  barplot(eig, names.arg=paste(round(eig, 2),"%"), las=2, main=eig_title)


  for (i in list.cpa){
    cpa_title=paste(cpa.title,": dim",i[1],"vs",i[2])
    plot.PCA(dat.PCA, choix=choix, axe=i, title=cpa_title, col.ind = col.ind)
    # second plot without label
    plot.PCA(dat.PCA, choix=choix, axe=i, title=cpa_title, col.ind = col.ind, label="none", lwd =4)
  }

  for (i in list.clust){
    if(dist.method == "euclidean"){
      dd <-  dist(dat.PCA[[choix]][["coord"]][,i], method = dist.method)
    }else{
      dd <- as.dist((1 - cor(t(dat.PCA[[choix]][["coord"]][,i]), method=dist.method))/2)
    }

    hc <- hclust(dd, method = "ward.D")
    clust_title <- paste("dendrogram on PCs projection :\nPCs",paste(i, collapse=", "), paste0("  dist : ",dist.method))

    mar.l <- max(nchar(hc$labels))/2


    hc <- as.dendrogram(hc)
    labels_colors(hc) <- col.ind[order.dendrogram(hc)]
    par(mar= c(mar.l,3,3,1))

    hc <- sort_hclust(hc)
    plot(hc, main = clust_title)
    # Second plot without label

    plot(hc ,labels=FALSE, main = clust_title, lwd = 4)
  }
  par(mar=old.par$mar)
}

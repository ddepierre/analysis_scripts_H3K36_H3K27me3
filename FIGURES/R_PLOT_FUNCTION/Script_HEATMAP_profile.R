######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################


# FUNCTION TO PLOT HEATMAP using matrix like output from DEEPTOOLS computeMatrix
## from https://github.com/al2na/genomation/blob/master/R/plotMatrix.R

######################################################################################################################################################
######################################################################################################################################################

######################################################################################################################################################
## LIBRARY
######################################################################################################################################################

library(grid)

######################################################################################################################################################
## FUNCTION
######################################################################################################################################################

.heatLegendY<-function(min,max,cols,legend.name,main=TRUE,cex.legend=1,
                       cex.lab=1){

  # get value range as a vector of 100
  vals=seq(min,max,length.out=100)
  rng <- range(vals, na.rm = TRUE) # get min/mqx
  m <- (vals - min)/(max-min) # get normalized range
  rasta= rgb(colorRamp(cols)(m), maxColorValue= 255) # get color for each element of range

  grid.raster( rev(rasta), interpolate=FALSE,height = unit(1, "npc"),
               width=unit(1, "npc")) # make the legend
  # make legend ticks
  at = seq(0,1,length.out=5); label = seq(min,max,length.out=5)

  #make the axis of the legend
  grid.yaxis(at=at,label=formatC(label,digits=2,format="g"),main=main,
             edits=gEdit("labels", rot=90,hjust=0.5),
             gp=gpar(cex=cex.legend)) # make axis for legend
  my.x = -2
  if(main==FALSE)
    my.x=3.4
  grid.text(legend.name,rot=90,x=unit(my.x, "npc"),gp=gpar(cex=cex.lab))

}

.convertToColors <- function(mat,cols, RangeValue = NULL) {
  # Produce 'normalized' version of matrix, with values ranging from 0 to 1
  # rng <- range(mat, na.rm = TRUE)
  if(is.null(RangeValue)){
    rng = range(mat, na.rm = TRUE)
    }else{
    rng = c(RangeValue[1], RangeValue[2])
  }
  m <- (mat - rng[1])/(diff(rng))
  # Convert to a matrix of sRGB color strings
  #m2 <- m; class(m2) <- "character"
  m2<-matrix("transparent",ncol=ncol(m),nrow=nrow(m))
  m2[!is.na(m)] <- rgb(colorRamp(cols)(m[!is.na(m)]), maxColorValue = 255)
  #m2[is.na(m)] <- "transparent"
  return(m2)
}

.gridHeat<-function(mat,col, RangeValue,xcoords,xlab,cex.lab,cex.axis,angle=0,
                    hjust=0,vjust=0){

  mat2=.convertToColors(mat,col, RangeValue)
  ras=grid.raster(mat2,interpolate = FALSE, width= unit(1, "npc"),
                  height=unit(1, "npc"))

  # make legend ticks
  at = seq(0,1,length.out=5); label = seq(min(xcoords),max(xcoords)
                                          ,length.out=5)

  ax=grid.xaxis(at=at,label=formatC(label,digits=4,format="g"),
                edits=gEdit("labels", rot=angle,hjust=hjust,vjust=vjust),
                gp=gpar(cex=cex.axis)) # make axis for legend

  grid.text(xlab,y=unit(-2.5, "lines"),gp=gpar(cex=cex.lab))
  #grid.draw(ax)
}

.rowSideCol<-function(group.vector,group.names=NULL,group.col=NULL,cex.lab=1){

  if( is.null(group.col) ){
    cols=rainbow(length(unique(group.vector)))
    img=cols[factor(group.vector,levels=unique(group.vector))]
  }else{
    img=group.col[factor(group.vector,levels=unique(group.vector))]
  }
  grid.raster(img,interpolate = FALSE, width= unit(1, "npc"),
              height=unit(1, "npc"))

  # segment heights calculated from group.vector
  # will be used to put group names in the middle of the segment
  segh=as.vector(table(factor(group.vector,levels=unique(group.vector))))
  name.coord=1-((cumsum(segh)-(segh/2))/sum(segh)) # NPC coord

  if( is.null(group.names)){
    grid.text(unique(group.vector), y=unit(name.coord,"npc"),
              x = unit(-0.5, "lines"),
              gp=gpar(cex=cex.lab),just="right")
  }else{
    grid.text(group.names, y=unit(name.coord,"npc"),
              x = unit(-0.5, "lines"),
              gp=gpar(cex=cex.lab),just="right")
  }

}

.heatLegendX<-function(min,max,cols,legend.name,main=TRUE,cex.legend=1,
                       cex.lab=1,hjust=0,vjust=0){

  # get value range as a vector of 100
  vals=seq(min,max,length.out=100)
  rng <- range(vals, na.rm = TRUE) # get min/mqx
  m <- (vals - min)/(max-min) # get normalized range
  rasta= rgb(colorRamp(cols)(m), maxColorValue = 255) # get color for each element of range

  grid.raster( matrix((rasta),nrow=1), interpolate=FALSE,height = unit(1, "npc"),
               width=unit(1, "npc")) # make the legend
  # make legend ticks
  label = pretty(c(min,max),n=5);at = seq(0,1,length.out=length(label));

  #make the axis of the legend
  grid.xaxis(at=at,label=label,main=main,
             edits=gEdit("labels",hjust=hjust,vjust=vjust),
             gp=gpar(cex=cex.legend)) # make axis for legend
  my.y = -3
  grid.text(legend.name,y=unit(my.y, "npc"),gp=gpar(cex=cex.lab))

}



.winsorize<-function(mat,rng){
  hi.th=quantile(mat,rng[2]/100,na.rm=TRUE)
  lo.th=quantile(mat,rng[1]/100,na.rm=TRUE)
  mat[mat>hi.th]=hi.th
  mat[mat<lo.th]=lo.th
  mat
}



heatMatrixMat=function (mat, grid = FALSE, col = NULL, xcoords = NULL, group = NULL,
    group.col = NULL, order = FALSE, RangeValue = NULL, winsorize = c(0, 100), kmeans = FALSE,
    k = 3, main = "", legend.name = NULL, cex.legend = 1, xlab = NULL,
    cex.main = 1, cex.lab = 1, cex.axis = 1, newpage = TRUE, colorSet = NULL)
{
require(gridBase)
if(is.null(colorSet)){
  .jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  # .jets<-colorRampPalette(c("#7F0000" ,   "red", "#FF7F00","yellow","#7FFF7F","cyan","#007FFF","blue", "#00007F"))
}else{
  .jets<-colorRampPalette(colorSet)
}
    if (class(mat) != "matrix") {
        stop("'mat' is not a matrix\n")
    }
    mat2 = mat
    group.vector = NULL
    group.names = NULL
    if(is.null(RangeValue)){
      if (winsorize[2] < 100 | winsorize[1] > 0) {
        hi.th = quantile(mat2, winsorize[2]/100, na.rm = TRUE)
        lo.th = quantile(mat2, winsorize[1]/100, na.rm = TRUE)
        mat2[mat2 > hi.th] = hi.th
        mat2[mat2 < lo.th] = lo.th
      }
    }else{

        mat2[mat2 > RangeValue[2]] = RangeValue[2]
        mat2[mat2 < RangeValue[1]] = RangeValue[1]
    }
    if (kmeans) {
        if (any(is.na(mat2))) {
            mat3 = impute.knn(mat2, k = 10, rowmax = 0.5, colmax = 0.8,
                maxp = 1500)$data
            clu = kmeans(mat3, c = k)
        }
        else {
            clu = kmeans(mat2, c = k)
        }
        group.vector = clu$cluster
        kcenters = clu$centers
        mat2 = mat2[order(group.vector), ]
        group.vector = group.vector[order(group.vector)]
        if (order) {
            g.factor = factor(group.vector, levels = unique(group.vector))
            cent.val = rowSums(kcenters, na.rm = TRUE)[g.factor]
            my.order = order(-cent.val, group.vector, -rowSums(mat2,
                na.rm = TRUE))
            mat2 = mat2[my.order, ]
            group.vector = group.vector[my.order]
        }
        group.names = unique(group.vector)
    }
    if (!is.null(group) & !kmeans) {
        if (is.list(group)) {
            win.numbs = (lapply(group, function(x) unique(x)))
            win.vec = unlist(win.numbs)
            if (any(table(unlist(win.numbs)) > 1)) {
                stop("'group' containing a list must not have duplicated numbers\n")
            }
            row.ids = rownames(mat2)
            group.vector = rep(0, length(row.ids))
            for (i in 1:length(win.numbs)) {
                group.vector[row.ids %in% win.numbs[[i]]] = i
            }
            if (!is.null(names(group))) {
                group.names = names(group)
            }
            if (all(group.vector == 0)) {
                stop("None of the elements in 'group' are a part of rownames(mat) \n")
            }
            if (any(group.vector == 0)) {
                warning("Number of elements in 'group' argument is less then nrow(mat) \n",
                  "Dropping rows from 'mat' that are not contained in 'group'\n")
                mat2 = mat2[group.vector > 0, ]
                group.vector = group.vector[group.vector > 0]
            }
        }
        else if (is.factor(group)) {
            if (length(group) != nrow(mat)) {
                stop("'group' is a factor, and its length should be equal to nrow(mat)\n")
            }
            group = factor(as.character(group), levels = as.character(unique(group)))
            group.names = levels(group)
            levels(group) = 1:length(levels(group))
            group.vector = as.numeric(group)
        }
        else {
            stop("'group' must be a factor or a list\n")
        }
        mat2 = mat2[order(group.vector), ]
        group.vector = group.vector[order(group.vector)]
        if (order) {
            my.order = order(group.vector, -rowSums(mat2, na.rm = TRUE))
            mat2 = mat2[my.order, ]
            group.vector = group.vector[my.order]
            names(group.vector) = rownames(mat2)
        }
    }
    else if (order & !kmeans) {
        order.vector = rep(1, nrow(mat2))
        names(order.vector) = rownames(mat2)
        order.vector = order.vector[order(-rowSums(mat2, na.rm = TRUE))]
        mat2 = mat2[order(-rowSums(mat2, na.rm = TRUE)), ]
    }
    if (!grid) {
        plot.new()
        vps <- baseViewports()
        pushViewport(vps$figure)
    }
    else {
        if (newpage)
            grid.newpage()
    }
    # X axis coordinates legend
    if (!is.null(xcoords) & is.vector(xcoords)) {
        if (length(xcoords) == 2 & xcoords[1] < xcoords[2]) {
            xcoords = seq(xcoords[1], xcoords[2], length.out = ncol(mat2))
        }
        if (length(xcoords) != ncol(mat2))
            stop("xcoords has wrong length: ", length(xcoords),
                " \n", " it should be equal to the number of columns of ScoreMatrix\n",
                " which is", ncol(mat2), "\n")
    }
    else {
        xcoords = 1:ncol(mat2)
    }
    if (is.null(col)) {
        col = .jets(100)
    }
    legendVp <- viewport(width = unit(0.7, "lines"), height = unit(0.4,
        "npc"), x = unit(0.71, "npc"), y = unit(0.5, "npc"),
        just = "left")
    pushViewport(legendVp)
    if(is.null(RangeValue)){
      rng = range(mat2, na.rm = TRUE)
      }else{
      rng = c(RangeValue[1], RangeValue[2])
    }
    .heatLegendY(min = rng[1], max = rng[2], col, legend.name,
        main = FALSE, cex.legend, cex.lab)
    popViewport()
    heatHeightNPC = 0.7
    heatVp <- viewport(width = unit(0.5, "npc"), height = unit(heatHeightNPC,
        "npc"), x = unit(0.2, "npc"), y = unit(0.5, "npc"), just = "left")
    pushViewport(heatVp)
    .gridHeat(mat2, col, RangeValue, xcoords, xlab, cex.lab, cex.axis, hjust = 0.5)
    popViewport()
    if (!is.null(group.vector)) {
        sideVp <- viewport(width = unit(0.05, "npc"), height = unit(heatHeightNPC,
            "npc"), x = unit(0.145, "npc"), y = unit(0.5, "npc"),
            just = "left")
        pushViewport(sideVp)
        grid.rect()
        .rowSideCol(group.vector, group.names = group.names,
            group.col = group.col, cex.lab = cex.lab)
        popViewport()
    }
    title.y = unit(0.9, "npc")
    grid.text(main, y = title.y, x = unit(0.45, "npc"), gp = gpar(cex = cex.main))
    if (!grid)
        popViewport()
    if (kmeans | !is.null(group.vector)) {
        return(invisible(group.vector))
    }
    else if (order & is.null(group.vector)) {
        return(invisible(order.vector))
    }
}

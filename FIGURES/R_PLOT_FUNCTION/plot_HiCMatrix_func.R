######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################


# Function that plot HiC 2D plot around a given region for 2 HiC matrices

# NOTE: this function is still under development, only paramters and code lines necessary to reproduce FIGURE 6 A and FIGURE SUPP 6 A are usable, others are commented

######################################################################################################################################################
######################################################################################################################################################


            library(gplots) # for textplot
            library("heatmap3")
            library(plyr)
            library(RColorBrewer)
            library(pheatmap)

            # plot_DiffHiCMatrix.func
                plot_DiffHiCMatrix.func = function(hic.mtx.WT, hic.mtx.KD, region = "2L:4000000-10000000", unzoom=NULL, res.nmb = res.nmb, main.fig = NULL, plot.name = NULL, plot.type = "svg", log2T=T, center.num = 0, colMinMax = c(0,10), center.num.mirror = 1, colMinMax.mirror = c(-1,1), outlier.cut=0.01, paletteLength.nmb = 50, palette.bias = 1, out.fig = out.fig){
                    # Helper
                        #Na. plot_HiCMatrix.func
                        #De. Function that plot HiC 2D plot around a given region
                        #Us. bedpe2Pairstbl.func(bedpe.gr = loop.bedpe.gr, res.nmb = res.nmb, apaSize.nmb = 41, minBinDist.nmb = 1, chromosomes.chr = chromosomes.chr, genome.chr="hg19")
                        #Ar. hic.mtx <"list" of matrices>: 
                        #Ar. region <char>: region to plot, eiher chr:start-end, either a GRanges object with 1 range (see exemple)
                        #Ar. res.nmb <integer>: The Hic Resolution
                        #Ar. main.fig <char>: title of the plot
                        #Ar. plot.name <char>: title of the pdf or svg
                        #Ar. out.fig <char>: out path of pdf or svg
                        #Va. A svg or pdf
                        #Ex. 
         

                    options(scipen=999)

                    if(plot.type %ni% "svg" & plot.type %ni% "pdf"){
                        stop("Error: plot.type must be one of 'svg' or 'pdf'.")
                    }

                    # Prepare HiC submatrix to plot
                    if(class(region) %in% "GRanges"){
                        chromToPlot = as.vector(seqnames(region))
                        chromMatrixToPlot = paste0(as.vector(seqnames(region)), "_",as.vector(seqnames(region)))
                        startAnchorToPlot =  start(region)
                        endAnchorToPlot = end(region)
                        startMatrixToPlot =  start(region)
                        endMatrixToPlot = end(region)
                        if(!is.null(unzoom)){
                            startMatrixToPlot =  startMatrixToPlot-unzoom
                            endMatrixToPlot =endMatrixToPlot+unzoom
                        }
                        if(is.null(plot.name)){
                            plot.name = paste0(main.fig, "_", chromToPlot, ":",startMatrixToPlot ,"-",endMatrixToPlot,"_",region$name)
                            cat("generating automatic plot.name as \n", plot.name,"\n") 
                        }
                    }else if(class(region) %in% "character"){
                        chromToPlot = unlist(strsplit(region, ":"))[1]
                        chromMatrixToPlot = paste0(unlist(strsplit(region, ":"))[1], "_",unlist(strsplit(region, ":"))[1])
                        startAnchorToPlot = as.numeric(unlist(strsplit(unlist(strsplit(region, ":"))[2], "-"))[1])
                        endAnchorToPlot = as.numeric(unlist(strsplit(unlist(strsplit(region, ":"))[2], "-"))[2])
                        startMatrixToPlot = as.numeric(unlist(strsplit(unlist(strsplit(region, ":"))[2], "-"))[1])
                        endMatrixToPlot = as.numeric(unlist(strsplit(unlist(strsplit(region, ":"))[2], "-"))[2])
                        if(!is.null(unzoom)){
                            startMatrixToPlot =  startMatrixToPlot-unzoom
                            endMatrixToPlot =endMatrixToPlot+unzoom
                        }
                        if(is.null(plot.name)){
                            plot.name = paste0(main.fig, "_", region)
                            cat("generating automatic plot.name as \n", plot.name,"\n") 
                        }
                    }

                    startMatrixToPlot.bin = floor(startMatrixToPlot/res.nmb)
                    endMatrixToPlot.bin = floor(endMatrixToPlot/res.nmb)
                    if(endMatrixToPlot<startMatrixToPlot){
                        stop("Error in region: start should be before end\nstart is ", startMatrixToPlot, "\n",
                                "end is ", endMatrixToPlot, "\n")
                    }
                    if(chromMatrixToPlot %ni% names(hic.mtx.WT)){
                        stop("Error in region: Submatrix name ", chromMatrixToPlot, " is not found in names(hic.mtx.WT)\n")
                    }

                    hic.WT.submtx = hic.mtx.WT[[chromMatrixToPlot]]
                    hic.WT.submtx = as.matrix(hic.WT.submtx[startMatrixToPlot.bin:endMatrixToPlot.bin,startMatrixToPlot.bin:endMatrixToPlot.bin])
                    hic.KD.submtx = hic.mtx.KD[[chromMatrixToPlot]]
                    hic.KD.submtx = as.matrix(hic.KD.submtx[startMatrixToPlot.bin:endMatrixToPlot.bin,startMatrixToPlot.bin:endMatrixToPlot.bin])

                    hic.mirror.submtx = hic.WT.submtx+t(hic.KD.submtx)

                    # hic.WT.submtx.smth = image.smooth(hic.WT.submtx)
                    # hic.KD.submtx.smth = image.smooth(hic.KD.submtx)
                    # hic.WT.submtx.smth = hic.WT.submtx.smth$z
                    # hic.KD.submtx.smth = hic.KD.submtx.smth$z

                    # hic.mirror.submtx.smth = hic.WT.submtx.smth+t(hic.KD.submtx.smth)
                    # colnames(hic.mirror.submtx.smth) = NULL
                    # rownames(hic.mirror.submtx.smth) = NULL

                    hic.WT.submtx[lower.tri(hic.WT.submtx, diag=F)] <- NA # set NA on lower triangle
                    hic.KD.submtx[lower.tri(hic.KD.submtx, diag=F)] <- NA # set NA on lower triangle
                    
                    


                    # if(log2T == T){

                    #     hic.submtx = log2(hic.KD.submtx+1)-log2(hic.WT.submtx+1)
                    # }else{
                    #     hic.submtx = hic.KD.submtx-hic.WT.submtx
                    # }

                    # hic.submtx.smth = image.smooth(hic.submtx)
                    # hic.submtx.smth = hic.submtx.smth$z

                    # ADD label anchors and chrom window
                    startAnchorToPlot.bin = floor(startAnchorToPlot/res.nmb)
                    endAnchorToPlot.bin = floor(endAnchorToPlot/res.nmb)
                    # colnames(hic.submtx) = rep("", dim(hic.submtx)[2])
                    # rownames(hic.submtx) = rep("", dim(hic.submtx)[1])
                    # colnames(hic.submtx)[dim(hic.submtx)[2]-(startAnchorToPlot.bin-startMatrixToPlot.bin)] = "*"
                    # rownames(hic.submtx)[endMatrixToPlot.bin-endAnchorToPlot.bin+1] = "*"

                    colnames(hic.mirror.submtx) = rep("", dim(hic.mirror.submtx)[2])
                    rownames(hic.mirror.submtx) = rep("", dim(hic.mirror.submtx)[1])
                    colnames(hic.mirror.submtx)[dim(hic.mirror.submtx)[2]-(startAnchorToPlot.bin-startMatrixToPlot.bin)] = "*"
                    rownames(hic.mirror.submtx)[endMatrixToPlot.bin-endAnchorToPlot.bin+1] = "*"

                    ########################################
                    ## DEF COLOR PALETTE
                    ########################################
                    red.diffHiC = colorRampPalette(colors=brewer.pal(6, 'Reds'), space='Lab',interpolate='spline',bias=palette.bias)
                    blue.diffHiC = colorRampPalette(colors=brewer.pal(6, 'Blues'), space='Lab',interpolate='spline',bias=palette.bias)

                    ########################################
                    ## Raw values centered + 
                    ########################################
                        # # Vectoriser la matrice pour le calcul des breaks
                        # vec = as.vector(hic.submtx)
                        # vec = vec[which(!is.na(vec))]
                        # # Breaks
                        #     brk.vec = compute_Breaks(xs = vec, quantile.lgk=F, n = paletteLength.nmb)
                        # # Colors
                        #     colors.vec = getColorPalette(brk.vec = brk.vec, 
                        #                                     palette.func.pos = red.diffHiC, 
                        #                                     palette.func.neg = blue.diffHiC, 
                        #                                     center.num = center.num,  center.col = "#ffffff", heatmap.func = "heatmap3")
                        #     length(brk.vec)
                        #     length(colors.vec)
                    ########################################
                    ## Raw values + remove outlier 1%
                    ########################################
                        # hic.submtx.rmOutlier = hic.submtx                        
                        # vec = as.vector(hic.submtx)
                        # vec = vec[which(!is.na(vec))]
                        # hic.submtx.rmOutlier[hic.submtx.rmOutlier < quantile(vec,outlier.cut)] <- quantile(vec,outlier.cut)
                        # hic.submtx.rmOutlier[hic.submtx.rmOutlier > quantile(vec,1-outlier.cut)] <- quantile(vec,1-outlier.cut)
                        # # Vectoriser la matrice pour le calcul des breaks
                        # vec.rmOutlier = as.vector(hic.submtx.rmOutlier)
                        # vec.rmOutlier = vec.rmOutlier[which(!is.na(vec.rmOutlier))]
                        # # Breaks
                        #     brk.vec.rmOutlier = compute_Breaks(xs = vec.rmOutlier, quantile.lgk=F, min.lim=quantile(vec,outlier.cut), max.lim=quantile(vec,1-outlier.cut), n = paletteLength.nmb)
                        # # Colors
                        #     colors.vec.rmOutlier = getColorPalette(brk.vec = brk.vec.rmOutlier, 
                        #                                     palette.func.pos = red.diffHiC, 
                        #                                     palette.func.neg = blue.diffHiC, 
                        #                                     center.num = center.num,  center.col = "#ffffff", heatmap.func = "heatmap3")
                        #     length(brk.vec.rmOutlier)
                        #     length(colors.vec.rmOutlier)
                    ########################################
                    ## Raw values + remove outlier 1% + centered on 1
                    ########################################
                        # hic.submtx.rmOutlier.centered = hic.submtx                        
                        # vec = as.vector(hic.submtx)
                        # vec = vec[which(!is.na(vec))]
                        # hic.submtx.rmOutlier.centered[hic.submtx.rmOutlier.centered < quantile(vec,outlier.cut)] <- quantile(vec,outlier.cut)
                        # hic.submtx.rmOutlier.centered[hic.submtx.rmOutlier.centered > quantile(vec,1-outlier.cut)] <- quantile(vec,1-outlier.cut)
                        # # Vectoriser la matrice pour le calcul des breaks
                        # vec.rmOutlier.centered = as.vector(hic.submtx.rmOutlier.centered)
                        # vec.rmOutlier.centered = vec.rmOutlier.centered[which(!is.na(vec.rmOutlier.centered))]
                        # # Breaks
                        #     brk.vec.rmOutlier.centered = compute_Breaks(xs = vec.rmOutlier.centered, quantile.lgk=F, center.num=center.num, min.lim=NULL, max.lim=NULL, n = paletteLength.nmb)
                        # # Colors
                        #     colors.vec.rmOutlier.centered = getColorPalette(brk.vec = brk.vec.rmOutlier.centered, 
                        #                                     palette.func.pos = red.diffHiC, 
                        #                                     palette.func.neg = blue.diffHiC, 
                        #                                     center.num = center.num,  center.col = "#ffffff", heatmap.func = "heatmap3")
                        #     length(brk.vec.rmOutlier.centered)
                        #     length(colors.vec.rmOutlier.centered)
                    ########################################
                    ## Raw values + remove colMinMax
                    ########################################
                        # hic.submtx.lim = hic.submtx
                        # hic.submtx.lim[hic.submtx.lim < colMinMax[1]] <- colMinMax[1]
                        # hic.submtx.lim[hic.submtx.lim > colMinMax[2]] <- colMinMax[2]
                        # # Vectoriser la matrice pour le calcul des breaks
                        # vec = as.vector(hic.submtx.lim)
                        # vec = vec[which(!is.na(vec))]
                        # # Breaks
                        #     brk.vec.lim = compute_Breaks(xs = vec, quantile.lgk=F, center.num=center.num, min.lim=colMinMax[1], max.lim=colMinMax[2], n = paletteLength.nmb)
                        # # Colors
                        #     colors.vec.lim = getColorPalette(brk.vec = brk.vec.lim, 
                        #                                     palette.func.pos = red.diffHiC, 
                        #                                     palette.func.neg = blue.diffHiC, 
                        #                                     center.num = center.num,  center.col = "#ffffff", heatmap.func = "heatmap3")
                        #     length(brk.vec.lim)
                        #     length(colors.vec.lim)

                    ########################################
                    ## Raw values + remove colMinMax mirror
                    ########################################
                        hic.mirror.submtx.lim = hic.mirror.submtx
                        hic.mirror.submtx.lim[hic.mirror.submtx.lim < colMinMax.mirror[1]] <- colMinMax.mirror[1]
                        hic.mirror.submtx.lim[hic.mirror.submtx.lim > colMinMax.mirror[2]] <- colMinMax.mirror[2]
                        # Vectoriser la matrice pour le calcul des breaks
                        vec = as.vector(hic.mirror.submtx.lim)
                        vec = vec[which(!is.na(vec))]
                        # Breaks
                            brk.vec.mirror.lim = compute_Breaks(xs = vec, quantile.lgk=F, center.num=center.num.mirror, min.lim=colMinMax.mirror[1], max.lim=colMinMax.mirror[2], n = paletteLength.nmb)
                        # Colors
                            colors.vec.mirror.lim = getColorPalette(brk.vec = brk.vec.mirror.lim, 
                                                            palette.func.pos = red.diffHiC, 
                                                            palette.func.neg = blue.diffHiC, 
                                                            center.num = center.num,  center.col = "#ffffff", heatmap.func = "heatmap3")
                            length(brk.vec.mirror.lim)
                            length(colors.vec.mirror.lim)
                    ########################################
                    ########################################
                    ## Raw values + remove colMinMax mirror
                    ########################################
                        # hic.mirror.submtx.smth.lim = hic.mirror.submtx.smth
                        # hic.mirror.submtx.smth.lim[hic.mirror.submtx.smth.lim < colMinMax.mirror[1]] <- colMinMax.mirror[1]
                        # hic.mirror.submtx.smth.lim[hic.mirror.submtx.smth.lim > colMinMax.mirror[2]] <- colMinMax.mirror[2]
                        # # Vectoriser la matrice pour le calcul des breaks
                        # vec = as.vector(hic.mirror.submtx.smth.lim)
                        # vec = vec[which(!is.na(vec))]
                        # # Breaks
                        #     brk.vec.mirror.smth.lim = compute_Breaks(xs = vec, quantile.lgk=F, center.num=center.num.mirror, min.lim=colMinMax.mirror[1], max.lim=colMinMax.mirror[2], n = paletteLength.nmb)
                        # # Colors
                        #     colors.vec.mirror.smth.lim = getColorPalette(brk.vec = brk.vec.mirror.lim, 
                        #                                     palette.func.pos = red.diffHiC, 
                        #                                     palette.func.neg = blue.diffHiC, 
                        #                                     center.num = center.num,  center.col = "#ffffff", heatmap.func = "heatmap3")
                        #     length(brk.vec.mirror.smth.lim)
                        #     length(colors.vec.mirror.smth.lim)
                    ########################################
                    ## PLOT
                    ########################################
                        if(plot.type %in% "svg"){
                            svg(file=paste0(out.fig,plot.name, ".svg"), bg="white")
                        }else if(plot.type %in% "pdf"){
                            pdf(paste0(out.fig,plot.name, ".pdf"), bg="white")
                        }

                            heatmap3(hic.mirror.submtx.lim,col=colors.vec.mirror.lim, breaks=brk.vec.mirror.lim, main = paste0("rm ColMinMax + cntrd"),
                                        Rowv=NA,Colv=NA,cex.main=0.5, lasCol = 2, xlab = paste0("anchors position : ", chromMatrixToPlot,":",startAnchorToPlot, "-",endAnchorToPlot, " (dist : ",endAnchorToPlot-startAnchorToPlot ,"bp)") , scale="none",margins = c(2, 2),revC=TRUE,
                                        balanceColor=F,legendfun=function() plot(0,0,bty="n",xaxt="n",yaxt="n",type="n"))

                            # heatmap3(hic.mirror.submtx.smth.lim,col=colors.vec.mirror.smth.lim, breaks=brk.vec.mirror.smth.lim, main = paste0("rm ColMinMax + cntrd"),
                            #             Rowv=NA,Colv=NA,cex.main=0.5, lasCol = 2, xlab = paste0("anchors position : ", chromMatrixToPlot,":",startAnchorToPlot, "-",endAnchorToPlot, " (dist : ",endAnchorToPlot-startAnchorToPlot ,"bp)") , scale="none",margins = c(2, 2),revC=TRUE,
                            #             balanceColor=F,legendfun=function() plot(0,0,bty="n",xaxt="n",yaxt="n",type="n"))

                            pheatmap(
                                t(as.matrix(brk.vec.mirror.lim)), 
                                cluster_rows=FALSE, cluster_cols=FALSE,
                                show_rownames=TRUE, show_colnames=TRUE,
                                color=colors.vec.mirror.lim,
                                breaks=brk.vec.mirror.lim,
                                fontsize=4, na_col="grey90", border_color=NA,
                                cellwidth = 1, cellheight = 100,
                                main=paste0(main.fig, " removed ColMinMax + centered / LEGEND COLOR")
                            )

                            # heatmap3(hic.submtx,col=colors.vec, breaks=brk.vec, main = "raw values",
                            #             Rowv=NA,Colv=NA,cex.main=0.5, lasCol = 2, xlab = paste0("anchors position : ", chromMatrixToPlot,":",startAnchorToPlot, "-",endAnchorToPlot, " (dist : ",endAnchorToPlot-startAnchorToPlot ,"bp)"), scale="none",margins = c(2, 2),revC=TRUE)

                            # heatmap3(hic.submtx.rmOutlier,col=colors.vec.rmOutlier, breaks=brk.vec.rmOutlier, main = paste0("rm ",outlier.cut*100,"% outliers"),
                            #             Rowv=NA,Colv=NA,cex.main=0.5, lasCol = 2, xlab = paste0("anchors position : ", chromMatrixToPlot,":",startAnchorToPlot, "-",endAnchorToPlot, " (dist : ",endAnchorToPlot-startAnchorToPlot ,"bp)"), scale="none",margins = c(2, 2),revC=TRUE)

                            # heatmap3(hic.submtx.rmOutlier.centered,col=colors.vec.rmOutlier.centered, breaks=brk.vec.rmOutlier.centered, main = paste0("rm ",outlier.cut*100,"% outliers + cntrd"),
                            #             Rowv=NA,Colv=NA,cex.main=0.5, lasCol = 2, xlab = paste0("anchors position : ", chromMatrixToPlot,":",startAnchorToPlot, "-",endAnchorToPlot, " (dist : ",endAnchorToPlot-startAnchorToPlot ,"bp)"), scale="none",margins = c(2, 2),revC=TRUE, 
                            #             balanceColor=F,legendfun=function() plot(0,0,bty="n",xaxt="n",yaxt="n",type="n"))
                            # pheatmap(
                            #     t(as.matrix(brk.vec.rmOutlier.centered)), 
                            #     cluster_rows=FALSE, cluster_cols=FALSE,
                            #     show_rownames=TRUE, show_colnames=TRUE,
                            #     color=colors.vec.rmOutlier.centered,
                            #     breaks=brk.vec.rmOutlier.centered,
                            #     fontsize=4, na_col="grey90", border_color=NA,
                            #     cellwidth = 1, cellheight = 100,
                            #     main=paste0(main.fig, " rm ",outlier.cut*100," % outliers + centered / LEGEND COLOR")
                            # )

                            # heatmap3(hic.submtx.lim,col=colors.vec.lim, breaks=brk.vec.lim, main = paste0("rm ColMinMax + cntrd"),
                            #             Rowv=NA,Colv=NA,cex.main=0.5, lasCol = 2, xlab = paste0("anchors position : ", chromMatrixToPlot,":",startAnchorToPlot, "-",endAnchorToPlot, " (dist : ",endAnchorToPlot-startAnchorToPlot ,"bp)") , scale="none",margins = c(2, 2),revC=TRUE,
                            #             balanceColor=F,legendfun=function() plot(0,0,bty="n",xaxt="n",yaxt="n",type="n"))

                            # heatmap3(hic.submtx.smth.lim,col=colors.vec.smth.lim, breaks=brk.vec.smth.lim, main = paste0("rm ColMinMax + cntrd"),
                            #             Rowv=NA,Colv=NA,cex.main=0.5, lasCol = 2, xlab = paste0("anchors position : ", chromMatrixToPlot,":",startAnchorToPlot, "-",endAnchorToPlot, " (dist : ",endAnchorToPlot-startAnchorToPlot ,"bp)") , scale="none",margins = c(2, 2),revC=TRUE,
                            #             balanceColor=F,legendfun=function() plot(0,0,bty="n",xaxt="n",yaxt="n",type="n"))

                            # pheatmap(
                            #     t(as.matrix(brk.vec.lim)), 
                            #     cluster_rows=FALSE, cluster_cols=FALSE,
                            #     show_rownames=TRUE, show_colnames=TRUE,
                            #     color=colors.vec.lim,
                            #     breaks=brk.vec.lim,
                            #     fontsize=4, na_col="grey90", border_color=NA,
                            #     cellwidth = 1, cellheight = 100,
                            #     main=paste0(main.fig, " removed ColMinMax + centered / LEGEND COLOR")
                            # )

                            # plot(density(c(hic.submtx), na.rm=T), main="contact density")



                        dev.off()
                        options(scipen=0)

                }
 




                # compute_Breaks
                    compute_Breaks <- function(xs, quantile.lgk=F, min.lim=NULL, max.lim=NULL, center.num=NULL, n = 10) {
                        # Helper
                            #Na. compute_Breaks
                            #De. Function that compute breaks from a vector for plot colors
                            #Us. compute_Breaks(xs = vec, quantile.lgk=F, min=NULL, max=NULL, n = paletteLength.nmb)
                            #Ar. xs vector of numeric values
                            #Ar. quantile.lgk quantile break or not
                            #Ar. min.lim min value for breaks
                            #Ar. max.lim max value for breaks
                            #Ar. center.num, if null, no centered, if set to a numeric, breaks will be centered on it. Only work without quantile
                            #Ar. n number of breaks
                            #Va. None
                            #Ex. vec = c(rnorm(100, 0, 1.5), rnorm(50, 4, 2))
                            #Ex. Bk = compute_Breaks(xs = vec, quantile.lgk=F, min=NULL, max=NULL, n = paletteLength.nmb)
                            #Ex. Bk
                            #Ex. Bk = compute_Breaks(xs = vec, quantile.lgk=F, min.lim=-3, max.lim=3, n = paletteLength.nmb)
                            #Ex. Bk
                            #Ex. Bk = compute_Breaks(xs = vec, quantile.lgk=T, min=NULL, max=NULL, n = paletteLength.nmb)
                            #Ex. Bk
                            #Ex. Bk = compute_Breaks(xs = vec, quantile.lgk=F, min=NULL, max=NULL, center.num = 0 , n = paletteLength.nmb)
                            #Ex. Bk
                            #Ex. 
                        # remove NA
                        xs = xs[which(!is.na(xs))]

                        if(quantile.lgk==T){
                            if(is.null(min.lim)){min.lim = min(xs)}
                            if(is.null(max.lim)){max.lim = max(xs)}
                            if(is.null(center.num)){
                                breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
                                breaks <- breaks[!duplicated(breaks)]
                            }else{
                                breaks.inf = quantile(xs[xs<center.num], probs = seq(0, 1, length.out = floor(n/2)))
                                breaks.sup = quantile(xs[xs>center.num], probs = seq(0, 1, length.out = floor(n/2)))
                                breaks <- c(breaks.inf, breaks.sup)
                                breaks <- breaks[!duplicated(breaks)]
                            }

                        }
                        if(quantile.lgk==F){
                            if(is.null(min.lim)){min.lim = min(xs)}
                            if(is.null(max.lim)){max.lim = max(xs)}
                            if(is.null(center.num)){
                                breaks <- quantile(c(min.lim, max.lim), probs = seq(0, 1, length.out = n))
                                breaks <- breaks[!duplicated(breaks)]
                            }else{
                                breaks <- c(seq(min.lim, center.num, length.out=ceiling(n/2) + 1),
                                             seq(center.num, max.lim, length.out=floor(n/2)))
                                breaks <- breaks[!duplicated(breaks)]
                            }
                        }
                        if(length(breaks)<5){cat("WARNING : breaks number <", length(breaks),"\n" )}
                        return(breaks)
                    }


                getColorPalette = function(brk.vec = brk.vec, palette.func.pos = red.diffHiC, palette.func.neg = blue.diffHiC, center.num = center.num, center.col="#ffffd8", heatmap.func = "heatmap3"){
                    if(heatmap.func == "heatmap3"){ # then length(break) == length(color)-1
                        if(is.null(center.num)){
                            if(length(which(brk.vec<0)) & length(which(brk.vec>0)) ){                                         # Si il y a des valeurs négatives ET positives
                                colors.vec = c(rev(palette.func.neg(length(which(brk.vec<0))-1)), center.col,palette.func.pos(length(which(brk.vec>0))-1))
                            }else if (length(which(brk.vec<0))){                                                              # Sinon, si il y a QUE des valeurs négatives
                                if(length(which(brk.vec<0)) == length(brk.vec)){
                                    colors.vec = rev(palette.func.neg(length(which(brk.vec<0))-1))
                                }else if(length(which(brk.vec<0)) == length(brk.vec)-1){
                                    colors.vec = rev(palette.func.neg(length(which(brk.vec<0))))
                                }
                            }else{  
                                if(length(which(brk.vec>0)) == length(brk.vec)){
                                    colors.vec = palette.func.pos(length(which(brk.vec>0))-1)
                                }else if(length(which(brk.vec>0)) == length(brk.vec)-1){
                                    colors.vec = palette.func.pos(length(which(brk.vec>0)))
                                }                                                                                   # Sinon, il y a QUE des valeurs positives
                            }
                        }else{
                            if(length(which(brk.vec<center.num)) & length(which(brk.vec>center.num)) ){                                         # Si il y a des valeurs négatives ET positives
                                colors.vec = c(rev(palette.func.neg(ceiling(length(brk.vec)/2)-1)), center.col,palette.func.pos(ceiling(length(brk.vec)/2)-1))
                            }else if (length(which(brk.vec<center.num))){                                                              # Sinon, si il y a QUE des valeurs négatives
                                if(length(which(brk.vec<center.num)) == length(brk.vec)){
                                    colors.vec = rev(palette.func.neg(length(which(brk.vec<center.num))-1))
                                }else if(length(which(brk.vec<center.num)) == length(brk.vec)-1){
                                    colors.vec = rev(palette.func.neg(length(which(brk.vec<center.num))))
                                }
                            }else{  
                                if(length(which(brk.vec>center.num)) == length(brk.vec)){
                                    colors.vec = palette.func.pos(length(which(brk.vec>center.num))-1)
                                }else if(length(which(brk.vec>center.num)) == length(brk.vec)-1){
                                    colors.vec = palette.func.pos(length(which(brk.vec>center.num)))
                                }                                                                                   # Sinon, il y a QUE des valeurs positives
                            }
                        }
                    }else{ # then length(break) == length(color) (as for pheatmap)
                        if(is.null(center.num)){
                            stop("getColorPalette not already coded")
                        }else{
                            if(length(which(brk.vec<center.num)) & length(which(brk.vec>center.num)) ){                                         # Si il y a des valeurs négatives ET positives
                                colors.vec = c(rev(palette.func.neg(ceiling(length(brk.vec)/2)-1)), center.col,palette.func.pos(ceiling(length(brk.vec)/2)))
                            }else if (length(which(brk.vec<center.num))){                                                              # Sinon, si il y a QUE des valeurs négatives
                                colors.vec = rev(palette.func.neg(length(which(brk.vec<center.num))))
                            }else{  
                                colors.vec = palette.func.pos(length(which(brk.vec>center.num)))
                            }
                        }
                    }
                    return(colors.vec)
                }

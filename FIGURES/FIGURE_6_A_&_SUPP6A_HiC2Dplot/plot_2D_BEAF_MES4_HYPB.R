######################################################################################################################################################
######################################################################################################################################################
# analysis_scripts_H3K36_H3K27me3
# Cuvier's Lab
# David Depierre
# dav.depierre@gmail.com
# 2022.07
######################################################################################################################################################
######################################################################################################################################################

workdir = "/home/depierre/Bureau/analysis_scripts_H3K36_H3K27me3/"

source(paste0(workdir, "FIGURES/R_PLOT_FUNCTION/plot_HiCMatrix_func.R"))



######################################################################################
######################################################################################
        ## LOAD HIC Matrix / here is only Chr2L HiC matrix, in an a list of sparse matrix format
                    HiC_luc_KD_rep1_inter_30_m_oe_n_VC_SQRT_r_1000 = readRDS(paste0(workdir, "DATA_PROCESSING/PROCESSED_HIC_FILES/HiC_luc_KD_rep1_inter_30_m_oe_n_VC_SQRT_r_1000_Chr2L.RDS"))
                    HiC_HYPB_KD_mergeRep_inter_30_m_oe_n_VC_SQRT_r_1000 = readRDS(paste0(workdir, "DATA_PROCESSING/PROCESSED_HIC_FILES/HiC_HYPB_KD_mergeRep_inter_30_m_oe_n_VC_SQRT_r_1000_Chr2L.RDS"))
                    HiC_MES4_KD_mergeRep_inter_30_m_oe_n_VC_SQRT_r_1000 = readRDS(paste0(workdir, "DATA_PROCESSING/PROCESSED_HIC_FILES/HiC_MES4_KD_mergeRep_inter_30_m_oe_n_VC_SQRT_r_1000_Chr2L.RDS"))


####################################################################################################
####################################################################################################

                plot_DiffHiCMatrix.func(hic.mtx.WT = HiC_luc_KD_rep1_inter_30_m_oe_n_VC_SQRT_r_1000,
                                        hic.mtx.KD = HiC_MES4_KD_mergeRep_inter_30_m_oe_n_VC_SQRT_r_1000,
                                        region = "2L:16650000-16750000", unzoom = NULL, res.nmb = 1000, 
                                        main.fig = "DIFFHiC_luc_KD_rep1_MES4_KD_mergeRep_inter_30_m_oe_n_VC_SQRT_r_1000",
                                        plot.name = NULL, plot.type = "pdf", log2T=F,
                                        center.num = 0, colMinMax = c(-10,10), center.num.mirror = 1, colMinMax.mirror = c(-8,3), outlier.cut=0.01,
                                        paletteLength.nmb = 50, palette.bias = 0.5,
                                        out.fig = paste0(workdir, "FIGURES/FIGURE_6_A_&_SUPP6A_HiC2Dplot/"))

                plot_DiffHiCMatrix.func(hic.mtx.WT = HiC_luc_KD_rep1_inter_30_m_oe_n_VC_SQRT_r_1000,
                                        hic.mtx.KD = HiC_HYPB_KD_mergeRep_inter_30_m_oe_n_VC_SQRT_r_1000,
                                        region = "2L:16650000-16750000", unzoom = NULL, res.nmb = 1000, 
                                        main.fig = "DIFFHiC_luc_KD_rep1_HYPB_KD_mergeRep_inter_30_m_oe_n_VC_SQRT_r_1000",
                                        plot.name = NULL, plot.type = "pdf", log2T=F,
                                        center.num = 0, colMinMax = c(-10,10), center.num.mirror = 1, colMinMax.mirror = c(-8,3), outlier.cut=0.01,
                                        paletteLength.nmb = 50, palette.bias = 0.5,
                                        out.fig = paste0(workdir, "FIGURES/FIGURE_6_A_&_SUPP6A_HiC2Dplot/"))


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
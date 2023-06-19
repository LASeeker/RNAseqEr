#' Quick and easy processing of scRNAseq data
#'
#' @description
#' This function takes a quality controlled, and if necessary batch corrected
#' and integrated Seurat object as input and performs a complete standard analysis
#' workflow including the finding of an appropriate clustering resolution,
#' differential gene expression analysis, and gene intology analysis.
#'
#' @param seur_obj quality controlled, and if necessary batch corrected
#' and integrated Seurat object
#' @param n_pcs number of dimensions (1:n_pcs) used for findNeighbors() and and
#' non linear dimensional reduction (UMAP, TSNE).
#' @param res resolutions for clustering
#' @param select_genes genes that should be scaled. Default is all genes in Seurat
#' object but setting to variable genes may be valid too.
#' @param elbow_dims How many dimensions should be plotted in elbow plot to
#' determine n_pcs. Default is
#' @param tsne TRUE/FALSE, default is true which means TSNE dim reduction is
#' performed additionally to UMAP. Default is FALSE.
#' @param col_pattern name pattern used to identify columns that contain clustering
#' information at different resolutions. Default is "RNA_snn_res."
#' @param plot_cols color palette default is colour_pallette() from the RNAseqEr
#' package
#' @param clust_lab Whether DimPlot should be annotated with clusters labels
#' @param label_size size of cluster labels
#' @param save_dir file path to where plots should be saved, default is working
#' directory
#' @param width plot width of pdf file for dimensional reduction plots. Default
#' is 7
#' @param height plot height of pdf file for dimensional reduction plots. Default
#' is 5
#' @param use_reduction which dimensional reduction should be used for plotting
#' ("pca". "tsne", "umap)
#' @param reduction_sil dimensional reduction that should be used for silhouette
#' plot. Default is "PCA"
#' @param clust_lab_sil TRUE/FALSE if should clisters be labeled in silhouette plot
#' @param width_sil width of silhouette plot, default is 7
#' @param height_sil height of output plot, default is 5
#' @param int_cols which colums should be used for differential gene expression
#' analysis. Default is all resolutions specified in res argument
#' @param only_pos whether only positive log-fold changes should be considered,
#' default is TRUE.
#' @param min_pct minimum proportion within cluster of interest expressing the
#' gene for differential gene expression analysis
#' @param logfc_threshold threshold for log2 fold changes
#' @param fil_pct_1 filtering threshold for the minimum proportion of cells within
#' the cluster of interest expressing the gene.
#' @param fil_pct_2 filtering threshold for the maximum proportion of cells outside
#' the cluster of interest expressing the gene.
#' @param test_use statistical approach being used for differential gene expression
#' analysis. Default is MAST.
#' @param int_cols_pw wich column(s) in the metadata containing clustering
#' information should be used for pairwise comparisons. More than one cluster
#' needs to be present at the chosen resoluton. Default resolution is "RNA_snn_res.1"
#' @param min_pct_pw minimum proportion within cluster of interest expressing the
#' gene for pairwise differential gene expression analysis. Default is 0.25
#' @param logfc_threshold_pw threshold for log2 fold changes used in pairwise
#' differential gene expression analysis. Default is 0.25
#' @param fil_pct_1_pw filtering threshold for the minimum proportion of cells within
#' the cluster of interest expressing the gene. Default is 0.25
#' @param fil_pct_2_pw iltering threshold for the maximum proportion of cells outside
#' the cluster of interest expressing the gene. Default is 0.1
#' @param assay_use Which assay should be used for DGE. Default is "RNA" to
#' ensure RNA assay is used for integrated objects.
#' @param ad_pval filtering parameter for compiling list of interesting marker genes.
#' Default filter for adjusted p value is 0.05 (smaller than).
#' @param avg_log filtering parameter for compiling list of interesting marker genes.
#' Default filter for adjusted p value is 1.2 (larger than).
#' @param pct_1 filtering parameter for compiling list of interesting marker genes.
#' Default filter for cells within the cluster of interest expressing the gene is 0.25.
#' @param pct_2 filtering parameter for compiling list of interesting marker genes.
#' Default filter for cells outside the cluster of interest expressing the gene is 0.1.
#' @param n_top number of top genes per cluster considered for compiling a list
#' of potentially interesting marker genes. Default is 10 genes per cluster.
#'
#' @return saves a lot of output data in newly generated folders. RETURNS AN
#' UPDATE SEURAT OBJ (IMPLEMENT!)
#' @export
#' @import Seurat SingleCellExperiment gtools ggsci RColorBrewer MAST
#'    here dplyr clusterProfiler org.Hs.eg.db SingleCellExperiment cluster
#'    bluster ggplot2 ggbeeswarm clustree
#' @examples
#' library(Seurat)
#' library(clustree)
#' seur <- quick_RNAseqEr(pbmc_small, tsne = FALSE)
quick_RNAseqEr <- function(seur_obj,
                           n_pcs = 20,
                           res = c(0.005, 0.01, 0.04, 0.05,
                                   seq(from = 0.1, to = 1, by = 0.1)),
                           select_genes = rownames(seur_obj),
                           elbow_dims = 50,
                           tsne = TRUE,

                           #for plotting dimensional reduction
                           col_pattern = "RNA_snn_res.",
                           plot_cols = colour_palette(),
                           clust_lab = TRUE,
                           label_size = 8,
                           save_dir = getwd(),
                           width=7,
                           height=5,
                           use_reduction = "umap",

                           # for silhouette plot
                           reduction_sil = "PCA",
                           clust_lab_sil = TRUE,
                           width_sil=7,
                           height_sil=5,

                           # perform differential gene expression at different
                           #clustering resolutions
                           int_cols = paste0(col_pattern, res),
                           only_pos = TRUE,
                           min_pct = 0.25,
                           logfc_threshold = 0.25,
                           fil_pct_1 = 0.25,
                           fil_pct_2 = 0.6,
                           test_use = "MAST",

                           # pairwise comparison of clusters at one or more
                           # resolutions of interest
                           int_cols_pw = "RNA_snn_res.1",
                           min_pct_pw = 0.25,
                           logfc_threshold_pw = 0.25,
                           fil_pct_1_pw = 0.25,
                           fil_pct_2_pw = 0.1,
                           assay_use = "RNA",

                           # compile list of interesting genes based on top x
                           # hits per cluster/group of interest using DGE
                           # results (pairwise and not)

                           ad_pval = 0.05,
                           avg_log = 1.2,
                           pct_1 = 0.25,
                           pct_2 = 0.6,
                           n_top = 10

                           ){
  # Perform standard Seurat processing
  seur_obj <- seurat_proc(seur_obj,
                          n_pcs = n_pcs,
                          res = res,
                          select_genes = select_genes,
                          elbow_dims = elbow_dims,
                          tsne = tsne)

  # plot dimensionally reduced data at different clustering resolutions and
  # save to file
  plot_dir_resol <- paste(save_dir,
                          "outs",
                          "plots",
                          "resolution_plots",
                          sep = "/")

  if(dir.exists(plot_dir_resol) == FALSE){
    dir.create(plot_dir_resol, recursive = TRUE)
    print("New directory created for saving plots in general and DimPlots at different
          resolutions in particular")
  }



  plot_list(seur_obj = seur_obj,
            col_pattern = col_pattern,
            plot_cols = plot_cols ,
            clust_lab = clust_lab,
            label_size = label_size,
            save_dir = plot_dir,
            width=width,
            height=height,
            use_reduction = use_reduction)

  # save cluster some indication of cluster stability
  sil_plot_dir <- paste(save_dir,
                        "outs",
                        "plots",
                        "silhouette_plot",
                        sep = "/")

  if(dir.exists(sil_plot_dir) == FALSE){
    dir.create(sil_plot_dir, recursive = TRUE)
    print("New directory created for saving silhouette plots")
  }


  sil_plot(seur_obj,
           reduction = reduction_sil,
           col_pattern = col_pattern,
           plot_cols = plot_cols,
           clust_lab = TRUE,
           save_to_file = TRUE,
           save_dir = sil_plot_dir,
           width = width_sil,
           height = height_sil)

  # approximate silhouette
  appr_sil_dir <- paste(save_dir,
                        "outs",
                        "plots",
                        "approx_silhouette",
                        sep = "/")

  if(dir.exists(appr_sil_dir) == FALSE){
    dir.create(appr_sil_dir, recursive = TRUE)
    print("New directory created for saving approximate silhouette plots")
  }


  approx_sil(seur_obj,
             reduction = reduction_sil,
             col_pattern = col_pattern,
             plot_cols = plot_cols,
             clust_lab = clust_lab_sil,
             label_size = label_size,
             save_dir = appr_sil_dir,
             height= height)

  # Calculate and plut cluster purity measures

  clu_pur_dir <- paste(save_dir,
                        "outs",
                        "plots",
                        "clu_pur_dir",
                        sep = "/")

  if(dir.exists(clu_pur_dir) == FALSE){
    dir.create(clu_pur_dir, recursive = TRUE)
    print("New directory created for saving cluster purity plots")
  }


  clu_pure(seur_obj,
           reduction = reduction_sil,
           col_pattern = col_pattern,
           plot_cols = plot_cols,
           clust_lab = clust_lab_sil,
           label_size = label_size,
           save_dir = clu_pur_dir,
           width=7,
           height=5)

  # Print cluster tree
  clu_tree_dir <- paste(save_dir,
                       "outs",
                       "plots",
                       "clu_tree_dir",
                       sep = "/")

  if(dir.exists(clu_tree_dir) == FALSE){
    dir.create(clu_tree_dir, recursive = TRUE)
    print("New directory created for saving a cluster tree plot")
  }

  cluster_tree <- clustree(seur_obj,
                           prefix = col_pattern,
                           exprs = c("data", "counts", "scale.data"),
                           assay = NULL
                           )

  pdf(paste(clu_tree_dir, "/cluster_tree.pdf"),
      paper="a4", width=8, height=11.5)

  print(cluster_tree)

  dev.off()



  #perform differential gene expression analysis at different resolutions and
  # save results to file
  dge_dir <- paste(save_dir,
                   "outs",
                   "tables",
                   "DGE",
                   "different_resol",
                   sep = "/")

  if(dir.exists(dge_dir) == FALSE){
    dir.create(dge_dir, recursive = TRUE)
    print("new directory created for saving output tables in general and results
          of differential gene expression analyses in particular")
  }


  all_res_mark <- int_res_all_mark(seur_obj = seur_obj,
                                   int_cols = int_cols,
                                   only_pos = only_pos,
                                   min_pct = min_pct,
                                   logfc_threshold = logfc_threshold,
                                   fil_pct_1 = fil_pct_1,
                                   fil_pct_2 = fil_pct_2,
                                   save_dir = dge_dir,
                                   test_use = test_use)

  # perform pairwise dge at chosen favourite resolution

  dge_pw_dir <- paste(save_dir,
                   "outs",
                   "tables",
                   "DGE",
                   "pairwise",
                   sep = "/")
  if(dir.exists(dge_pw_dir) == FALSE){
    dir.create(dge_pw_dir, recursive = TRUE)
    print("new directory created for saving results
          of pariwise differential gene expression analyses in particular")
  }

  pw_mark <- pairwise_dge(seur_obj = seur_obj,
                           int_cols = int_cols_pw,
                           only_pos = only_pos,
                           min_pct = min_pct_pw,
                           logfc_threshold = logfc_threshold_pw,
                           fil_pct_1 = fil_pct_1_pw,
                           fil_pct_2 = fil_pct_2_pw,
                           save_dir = dge_pw_dir,
                           test_use = "MAST",
                           assay_use = "RNA")

  # use DGE results to plot heatmap and decide on clustering resolution


  clu_mark <- gen_mark_list(file_dir = dge_dir,
                            ad_pval = ad_pval,
                            avg_log = avg_log,
                            pct_1 = pct_1,
                            pct_2 = pct_2,
                            pairwise = FALSE,
                            n_top = n_top)

  clu_mark_pw <- gen_mark_list(file_dir = dge_pw_dir,
                            ad_pval = ad_pval,
                            avg_log = avg_log,
                            pct_1 = pct_1,
                            pct_2 = pct_2,
                            pairwise = TRUE,
                            n_top = n_top)

  # concatenate genes of clu_mark and clu_mark_pw to list of interesting genes

  int_genes <- c(clu_mark$gene, clu_mark_pw$gene)
  int_genes <- unique(int_genes)


  #plot heatmap with interesting genes





}

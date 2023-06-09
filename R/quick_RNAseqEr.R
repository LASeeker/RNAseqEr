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
#' @param use_resol For heatmap. TRUE/FALSE whether column labels should be used that indicate
#' standard names for different clustering resolutions opposed to user indicated
#' cluster labels that for example specify a cell type. Default is TRUE.
#' @param col_pattern_hm For heatmap. If use_resol = TRUE, the column name pattern
#' that should be used for heatmap. Likely to equal col_pattern above.
#' @param col_names For heatmap. If use_resol = FALSE a string of at least one column names
#' has to be provided.
#' @param label_hm TRUE/FALSE whether cluster labels should be displayed on heatmap.
#' Default is TRUE
#' @param draw_lines For heatmap. TRUE/FALSE whether white lines should be drawn between columns
#' corresponding to each cluster. Default is FALSE
#' @param min_log2FC_go min log2FC considered for GO. Default is 0.25.
#' @param reverse_go TRUE/FALSE whether ADD DESCRIPTION. Default is FALSE.
#' @param translate_gene_id_from format in which genes are saved in Seurat object,
#' default is "SYMBOL".
#' @param translate_gene_id_to Format to which gebe symbol should be translated
#' for gene ontology analysis. Default is "ENTREZID"
#' @param org_use organism being used for gene ontology. Default is "org.Hs.eg.db"
#' which is human
#' @param ontology What kind of gene ontology should be performed. Default is "BP" for
#' biological pathway but "MF" for molecular function for example is valid, too.
#' @param pvalue_cutoff_go P-value cutoff for gene ontology. Default is 0.05.
#' @param qvalue_cutoff_go Cutoff for gene ontology q-value. Default is 0.05.
#' @param read_able TRUE/FALSE whether gene names should be readable in output tables
#' SAVE OUTPUT TABLES AND PLOTS!
#' @param do_seurat_proc TRUE/FALSE (default = TRUE) whether standard Seurat processing including
#' normalisation, dimensional reduction and clustering at different resolutions
#' should be performed.
#' @param do_dimplots TRUE/FALSE (default = TRUE) whether DimPlots at different resolutions should
#' be plotted and saved to file.
#' @param do_silhouette TRUE/FALSE (default = TRUE) whether silhouette plots
#' should be plotted and saved to file
#' @param do_purity TRUE/FALSE (default = TRUE) whether cluster purity plots
#' should be plotted and saved to file
#' @param do_dge TRUE/FALSE (default = TRUE) whether differential gene
#' expression should be performed. N.B. this may take a long time as
#' DGE will be performed at all resolutions and also at pairwise for all possible
#' cluster combinations at all resolutions
#' @param do_heatmap TRUE/FALSE (default = TRUE) whether heatmaps
#' should be plotted and saved to file
#' @param do_go TRUE/FALSE (default = TRUE) whether gene ontology should be performed
#' THINK ABOUT GENE INPUT AND GROUPING
#'
#' @return saves output data in newly generated folders. Returns an
#' updated Seurat object with clustering information and dimensional reductions.
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
                           n_top = 10,

                           #for creating heatmaps
                           use_resol = TRUE,
                           col_names,
                           col_pattern_hm = "RNA_snn_res",
                           label_hm = TRUE,
                           draw_lines = FALSE,

                           #for running gene ontology analysis
                           gene_list,
                           min_log2FC_go = 0.25,
                           reverse_go = FALSE,
                           translate_gene_id_from = "SYMBOL",
                           translate_gene_id_to = "ENTREZID",
                           org_use = "org.Hs.eg.db",
                           ontology = "BP",
                           pvalue_cutoff_go = 0.05,
                           qvalue_cutoff_go = 0.05,
                           read_able = TRUE,


                           #specify which steps to do
                           do_seurat_proc = TRUE,
                           do_dimplots = TRUE,
                           do_silhouette = TRUE,
                           do_purity = TRUE,
                           do_dge = TRUE,
                           do_heatmap = TRUE,
                           do_go =FALSE #SET TO TRUE IF CORRECTLY IMPLEMENTED
                           ){
  # Perform standard Seurat processing
  if(do_seurat_proc == TRUE){
  seur_obj <- seurat_proc(seur_obj,
                          n_pcs = n_pcs,
                          res = res,
                          select_genes = select_genes,
                          elbow_dims = elbow_dims,
                          tsne = tsne)
  }

  # plot dimensionally reduced data at different clustering resolutions and
  # save to file
  if(do_dimplots == TRUE){
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
            save_dir = plot_dir_resol,
            width=width,
            height=height,
            use_reduction = use_reduction)
  }

  # save cluster some indication of cluster stability
  if(do_silhouette == TRUE) {
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
  }

  # Calculate and plot cluster purity measures
  if(do_purity == TRUE){

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
  }

  #read in purity measures to see which cluster resolution is the one with the
  #largest number of clusters and largest purity measure
  setwd(clu_pur_dir)
  pur_dat_dir <- paste("../../../outs",
                       "tables",
                       "cluster_purity_data",
                       sep = "/")

  files <- list.files(pur_dat_dir, pattern = ".csv")
  myfiles <- lapply(paste(pur_dat_dir, files, sep = "/"), read.csv)

  sum_pure_df <- data.frame()

  for (j in 1:length(myfiles)) {
    dat <- myfiles[[j]]
    clu_count <- length(levels(as.factor(dat$cluster)))
    clu_pur_mean <- mean(dat$purity)
    pure_measure <- (clu_pur_mean^22)/clu_count
    df_temp <- data.frame(clu_res <- files[j],
                          mean_pur = clu_pur_mean,
                          clu_count = clu_count,
                          pure_measure =  pure_measure)


    if(j == 1){
      keep_df <- df_temp
      keep_res <- files[j]
      keep_num_clu<- clu_count
    }else{
      keep_df <- rbind(keep_df, df_temp)
      if(clu_pur_mean > 0.99 & clu_count >keep_num_clu){
        keep_res <- files[j]
        keep_num_clu<- clu_count
      }
    }
  }

    keep_df$keep_res <- rep(keep_res, nrow(keep_df))
    keep_df$keep_val <- rep(keep_val, nrow(keep_df))


    keep_resolution <- keep_df$keep_res[1]
    keep_resolution <- strsplit(keep_resolution, "_clu")[[1]][1]

    Idents(seurat_obj) <- keep_resolution



  write.csv(keep_df, paste0(pur_dat_dir, "/summary_purity_stats.csv"))

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

  pdf(paste0(clu_tree_dir, "/cluster_tree.pdf"),
      paper="a4", width=8, height=11.5)

  print(cluster_tree)

  dev.off()



  #perform differential gene expression analysis at different resolutions and
  # save results to file
  if(do_dge == TRUE){
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
  }


  #plot heatmap with interesting genes
  if(do_heatmap == TRUE){
  hm_dir <- paste(save_dir,
                  "outs",
                  "plots",
                  "heatmaps",
                  sep = "/")

  if(dir.exists(hm_dir) == FALSE){
    dir.create(hm_dir, recursive = TRUE)
    print("New directory created for saving heatmaps")
  }

  heatmap_seqEr(seur_obj,
                use_resol = use_resol,
                col_names,
                col_pattern = col_pattern_hm,
                int_genes,
                save_dir = hm_dir,
                label = label_hm,
                plot_cols = colour_palette(),
                draw_lines = draw_lines)
  }

  if(do_go == TRUE){
    perform_go(seur_obj,
               gene_list,
               min_log2FC = min_log2FC_go,
               reverse = reverse_go,
               translate_gene_id_from = translate_gene_id_from,
               translate_gene_id_to = translate_gene_id_to,
               org_use = org_use,
               ontology = ontology,
               pvalue_cutoff = pvalue_cutoff_go,
               qvalue_cutoff = qvalue_cutoff_go,
               read_able = read_able)

  }

  return(seur_obj)



}

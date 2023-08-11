#' Differential Abundance testing using Milo
#' @description
#' This function performs an analysis to see if there are more or less cells/ nuclei
#' per cluster depending on a metadata factor such as condition.
#' @param seur_obj quality controlled and clustered seurat object
#' @param cluster_col name of metadata column that should be considered for
#' clustering. For example "rough_annot" in the cns dataset from the RNAseqEr
#' library.
#' @param sample_id unique sample id which is a metadata column in seur_obj that
#' contains sample ids such as "uniq_id" in the RNAseqEr cns dataset.
#' @param cols_interest metadata columns the user would like to check. Those can be
#' factors of interest such as disease vs. control or young vs. old and can
#' should also include sample id and donor id for example. They shouls also include
#' all columns that are to be used to correct models by being used in account_for.
#' Examples may be c("uniq_id", "Tissue", "gender", "AgeGroup", "caseNO",
#' "SequencingPool", "X10XBatch") in the
#' RNAseqEer cns dataset.
#' @param test_factors  metadata columns that should be used for differential
#' abundance. The easiest test is between two groups such as disease and control.
#' If there are more than 2 levels, a pairwise comparison will be run (therefore,
#' the selection of factors with many levels is not ideal here.)
#' An example is c("Tissue", "gender", "AgeGroup") in the RNAseqEer cns dataset.
#' @param account_for add factors to statistically correct for in the model if required.
#' The default is FALSE where no additional factors are acconted for.
#' If model is to be corrected for factors, they have to be column names in the
#' seur_obj metadata and appear in the design matrix
#' (they have to be have to be element of cols_interest but not of test_factors).
#' They need to be supplied as strings concatenated and ending with a + symbol like this:
#' "X10XBatch +" or "X10XBatch + SequencingPool +").
#' @param k An integer scalar that specifies the number of nearest-neighbours to consider for the graph building.
#' @param d he number of dimensions to use if the input is a matrix of cells X
#' reduced dimensions. If this is provided, transposed should also be set=TRUE.
#' @param reduced_dim Dimensional reduction used for Milo. Default is "PCA"
#' @param refined 	A logical scalar that determines the sampling behavior used
#' by Milo. Default=TRUE implements a refined sampling scheme, specified by the
#' refinement_scheme argument.
#' @param prop Milo: A double scalar that defines what proportion of graph
#' vertices to randomly sample. Must be 0 < prop < 1. Default is 0.2.
#' @param plot_dimred Dimensional reduction used for plot. Default is "UMAP"
#' @param plot_text size of plot text. Default is 3.
#' @param plot_point size for points in plot. Default is 0.5
#' @param plot_alpha Milo: significance level for Spatial FDR (default: 0.1)
#' @param save_dir directory that is to be used as root to save output files.
#' Default is getwd()
#' @param dir_lab Celltype label that specifies if all cell types are
#' being tested (Default = "all_celltypes") or a cell lineage which is important
#' for structuring output in a systematic way.
#' @param height Height used for output plots in pdf files. Default is 5
#' @param width Width used for output plots in pdf files. Default is 8
#' @param beeswarm_colour vector of three colours that are used for beeswarm
#' plots. Default is c("red", "blue", "grey").
#'
#' @return good question
#' @import Seurat SingleCellExperiment miloR ggplot2 statmod scater patchwork here
#' @export
#'
#' @examples
#' milo_obj <- abundance_test(seur_obj = cns,
#'                       cluster_col = "rough_annot",
#'                       sample_id = "uniq_id",
#'                       cols_interest = c("uniq_id",
#'                                          "Tissue",
#'                                          "gender",
#'                                          "AgeGroup",
#'                                          "caseNO"),
#'                       test_factors = c("Tissue", "gender", "AgeGroup"),
#'                       save_dir = getwd()
#'

abundance_test <- function(seur_obj,
                           cluster_col,
                           sample_id,
                           cols_interest,
                           test_factors,
                           account_for = FALSE, # use for example "AgeGroup + " # has to be in cols_interest
                           k = 30,
                           d = 30,
                           reduced_dim = "PCA",
                           refined = TRUE,
                           prop = 0.2,
                           plot_dimred = "UMAP",
                           plot_text = 3,
                           plot_point = 0.5,
                           plot_alpha = 0.1,
                           save_dir = getwd(),
                           dir_lab = "all_celltypes",
                           height = 5,
                           width = 8,
                           beeswarm_colour = c("red", "blue", "grey")
                           ){
  milo_meta <- seur_obj@meta.data
  milo_obj <- Milo(as.SingleCellExperiment(seur_obj))

  milo_obj <- buildGraph(milo_obj,
                         k = k,
                         d = d,
                         reduced.dim = reduced_dim)

  milo_obj <- makeNhoods(milo_obj, k = k, d = d, refined = refined, prop = prop,
                         reduced_dims = reduced_dim)

  milo_plot <- plotNhoodSizeHist(milo_obj) +
    geom_vline(xintercept = 50, linetype="dotted",
               color = "blue", linewidth=1.5) +
    geom_vline(xintercept = 100, linetype="dotted",
               color = "red", linewidth=1.5)

  print(milo_plot)

  milo_dir_plot <- paste0(save_dir,
                          "/outs/",
                          dir_lab,
                          "/plots/miloR")
  if(dir.exists(milo_dir_plot) == FALSE){
    dir.create(milo_dir_plot, recursive = TRUE)
  }


  milo_dir_tab <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/tables/miloR")
  if(dir.exists(milo_dir_tab) == FALSE){
    dir.create(milo_dir_tab, recursive = TRUE)
  }


  milo_dir_dat <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/data/miloR")
  if(dir.exists(milo_dir_dat) == FALSE){
    dir.create(milo_dir_dat, recursive = TRUE)
  }

  pdf(paste0(milo_dir_plot,
             "/neighbourhood_size_histogram.pdf"), width=width, height=height)
  plot(milo_plot)

  dev.off()


  milo_obj <- countCells(milo_obj, samples = sample_id, meta.data = milo_meta)

  milo_design <- data.frame(colData(milo_obj)[, cols_interest])



  milo_design <- distinct(milo_design)
  rownames(milo_design) <- milo_design[[sample_id]] # that does not work. needed??

  milo_obj <- calcNhoodDistance(milo_obj, d = d, reduced.dim = reduced_dim)

  write.csv(milo_design, paste0(milo_dir_tab, "milo_design_df.csv"))


  #Testing
  for(i in 1: length(test_factors)){
    curr_fact <- test_factors[i]

    if(length(levels(as.factor(milo_meta[[curr_fact]]))) == 2){

      if(account_for != FALSE){
        da_results <- testNhoods(milo_obj, design = as.formula(paste0("~ ",
                                                                      account_for,
                                                                      curr_fact)),
                                 design.df = milo_design)

      }else{

      da_results <- testNhoods(milo_obj, design = as.formula(paste0("~ ",
                                                                    curr_fact)),
                               design.df = milo_design)
      }

      print(head(da_results))
      da_results %>%
        arrange(SpatialFDR) %>%
        head()

      qc_plot <- ggplot(da_results, aes(PValue)) + geom_histogram(bins=50) +
        ggtitle(curr_fact)

      print(qc_plot)

      pdf(paste0(milo_dir_plot,
                 "/",
                 curr_fact,
                 "_P-val_histogram.pdf"), width=width, height=height)
      plot(qc_plot)

      dev.off()

      qc_plot2 <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
        geom_point() +
        geom_hline(yintercept = 1)  +
        ggtitle(curr_fact)

      print(qc_plot2)

      pdf(paste0(milo_dir_plot,
                 "/",
                 curr_fact,
                 "_FDR_logFC.pdf"), width=width, height=height)
      plot(qc_plot)

      dev.off()

      milo_obj <- buildNhoodGraph(milo_obj)

      umap_pl <- plotReducedDim(milo_obj,
                                dimred = plot_dimred,
                                colour_by= curr_fact,
                                text_by = cluster_col,
                                text_size = plot_text, point_size = plot_point) +
        guides(fill="none")  +
        ggtitle(curr_fact)

      nh_graph_pl <- plotNhoodGraphDA(milo_obj,
                                      da_results,
                                      layout = plot_dimred,
                                      alpha = plot_alpha)


      save_plot <- umap_pl + nh_graph_pl +
        plot_layout(guides="collect")

      print(save_plot)

      pdf(paste0(milo_dir_plot,
                 "/",
                 curr_fact,
                 "_umap.pdf"), width=width*2, height=height)
      plot(save_plot)

      dev.off()


      da_results <- annotateNhoods(milo_obj,
                                   da_results,
                                   coldata_col = cluster_col)
      head(da_results)


      da_results$clust <- ifelse(da_results[paste0(cluster_col, "_fraction")] < 0.7,
                                 paste("Mixed"),
                                 paste(da_results[[paste0(cluster_col, "_fraction")]]))

      da_results$fdr <- da_results$FDR

      da_results$FDR <- ifelse(da_results$fdr < 0.05, "0.05",
                               ifelse(da_results$fdr < 0.1, "0.1", "> 0.1"))

      da_results$FDR <- factor(da_results$FDR, levels = c("0.05", "0.1",  "> 0.1"))

      beeswarm <- ggplot(da_results, aes(x = .data[[cluster_col]], y = logFC , color = FDR)) +
        geom_quasirandom()+
        scale_color_manual(values = beeswarm_colour) + coord_flip()+
        theme_minimal() +
        xlab("Cluster ID") +
        ylab(paste0("(<- ",
                    levels(as.factor(seur_obj@meta.data[[curr_fact]]))[1],
                    ") Log Fold Change (",
                    levels(as.factor(seur_obj@meta.data[[curr_fact]]))[2],
                    " -->)")) +
        ggtitle(curr_fact)


      pdf(paste0(milo_dir_plot,
                 "/",
                 curr_fact,
                 "_",
                 "_beeswarm.pdf"), width=width*2, height=height)
      plot(beeswarm)

      dev.off()



      plot(beeswarm)



      write.csv(da_results, paste0(milo_dir_tab, "/diff_abund_", curr_fact, ".csv"))


    } else if(length(levels(as.factor(milo_meta[[curr_fact]]))) > 2){
      # do all pairwise comparisons
      all_levels <- levels(as.factor(milo_meta[[curr_fact]]))
      all_comb <- combn(all_levels, 2)


      for(k in 1: ncol(all_comb)){
        da_results <- testNhoods(milo_obj,
                                 design = as.formula(paste0("~ 0 +",
                                                   curr_fact)),
                                 design.df = milo_design,
                                 model.contrasts = paste0(curr_fact,
                                                          all_comb[2,k],
                                                          " - ",
                                                          curr_fact,
                                                          all_comb[1,k]))


        da_results %>%
          arrange(SpatialFDR) %>%
          head()


        ggplot(da_results, aes(PValue)) +
          geom_histogram(bins=50) +
          ggtitle(paste0(curr_fact,
                         " ",
                         all_comb[2,k],
                         " vs. ",
                         all_comb[1,k]))

        ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
          geom_point() +
          geom_hline(yintercept = 1) +
          ggtitle(paste0(curr_fact,
                         " ",
                         all_comb[2,k],
                         " vs. ",
                         all_comb[1,k]))

        milo_obj <- buildNhoodGraph(milo_obj)
        umap_pl <- plotReducedDim(milo_obj,
                                  dimred = plot_dimred,
                                  colour_by= curr_fact,
                                  text_by = cluster_col,
                                  text_size = plot_text, point_size = plot_point) +
          guides(fill="none") +
          ggtitle(paste0(curr_fact,
                         " ",
                         all_comb[2,k],
                         " vs. ",
                         all_comb[1,k]))


        nh_graph_pl <- plotNhoodGraphDA(milo_obj,
                                        da_results,
                                        layout = plot_dimred,
                                        alpha = plot_alpha)

        save_plot <- umap_pl + nh_graph_pl +
          plot_layout(guides="collect")
        print(save_plot)


        pdf(paste0(milo_dir_plot,
                   "/",
                   curr_fact,
                   "_umap.pdf"), width=width*2, height=height)
        plot(save_plot)

        dev.off()

        da_results <- annotateNhoods(milo_obj,
                                     da_results,
                                     coldata_col = cluster_col)
        print(head(da_results))

        da_results$clust <- ifelse(da_results[paste0(cluster_col, "_fraction")] < 0.7,
                                   paste("Mixed"),
                                   paste(da_results[[paste0(cluster_col, "_fraction")]]))

        da_results$fdr <- da_results$FDR

        da_results$FDR <- ifelse(da_results$fdr < 0.05, "0.05",
                                 ifelse(da_results$fdr < 0.1, "0.1", "> 0.1"))

        da_results$FDR <- factor(da_results$FDR, levels = c("0.05", "0.1",  "> 0.1"))

        beeswarm <- ggplot(da_results, aes(x = .data[[cluster_col]], y = logFC , color = FDR)) +
          geom_quasirandom()+
          scale_color_manual(values = beeswarm_colour) + coord_flip()+
          theme_minimal() +
          xlab("Cluster ID") +
          ylab(paste0("(<- ",
                      all_comb[1,k],
                      ") Log Fold Change (",
                      all_comb[2,k],
                      " -->)")) +
          ggtitle(curr_fact)

        plot(beeswarm)

        pdf(paste0(milo_dir_plot,
                   "/",
                   curr_fact,
                   "_",
                   all_comb[2,k],
                   "_vs_",
                   all_comb[1,k],
                   "_beeswarm.pdf"), width=width*2, height=height)
        plot(beeswarm)

        dev.off()


        write.csv(da_results, paste0(milo_dir_tab,
                                     "/diff_abund_",
                                     curr_fact,
                                     "_",
                                     all_comb[2,k],
                                     "_vs_",
                                     all_comb[1,k], ".csv"))



      }


    } else {
      print("Choose factor with more than 1 level to test.")
    }


    saveRDS(milo_obj, paste0(milo_dir_dat,
                             "/milo_object_",
                             dir_lab,
                             ".RDS"))


  }
  return(milo_obj)

}










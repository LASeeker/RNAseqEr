#' Save feature plots
#' @description
#' Save feature plots generated from gene lists to file
#'
#'
#' @param seur_obj seurat object
#' @param marker_list list of marker genes that should be plotted
#' @param dir_lab label that indicated cell lineage to structure output. The
#' default is "all_celltypes"
#' @param save_label label that is provided to structure output and clarifies
#' if condition markers or cluster markers are being plotted.
#' @param save_dir directory where the output should be saved. Default is the
#' current working directory.
#' @param numb_genes number of genes that should be plotted together. Default is 16.
#' @param plotheight height of plot.Default is 18.
#' @param plotwidth width of plot. default is 20.
#' @param split_by If feature plot should be plotted separately for a condition
#' for example the corresponding metadata column name can be provided here.
#' Otherwise the default is "NULL" and the plots will not be splitted.
#' @param n_col number of columns used in output plot grid. The default is 4.
#'
#' @return Saves plot to file.
#' @import Seurat
#' @export
#'
#' @examples
#' markers <- c("SNAP25", "PLP1", "MAG", "CD74")
#' save_feat_plots(cns, marker_list = markers)
save_feat_plots <- function(seur_obj,
                            marker_list,
                            dir_lab = "all_celltypes",
                            save_label = "cluster_mark",
                            save_dir = getwd(),
                            numb_genes = 9,
                            plotheight = 18,
                            plotwidth = 20,
                            split_by = "NULL",
                            n_col = 4){
  gene_groups <- split(marker_list,
                       ceiling(seq_along(marker_list) / numb_genes))

  for(i in 1:length(gene_groups)){
    if(split_by != "NULL"){
      feature_plot<-FeaturePlot(seur_obj,
                                features = unlist(gene_groups[i]),
                                ncol = n_col, split.by = split_by)
    }else{
      feature_plot<-FeaturePlot(seur_obj,
                                features = unlist(gene_groups[i]),
                                ncol = n_col)
    }

    dir_save <- paste0(save_dir,
                       "/outs/",
                       dir_lab,
                       "/plots/",
                       save_label,
                       "/FeaturePlots")
    if(dir.exists(dir_save) == FALSE){
      dir.create(dir_save, recursive = TRUE)
    }


    pdf(paste0(dir_save, "/", save_label, "_", i, ".pdf", sep=""),
        height=plotheight, width = plotwidth)

    print(feature_plot)

    dev.off()

    print(feature_plot)

  }
}

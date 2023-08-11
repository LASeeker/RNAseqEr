#' Plotting violin plots from gene lists
#' @description
#' A gene list that should be filtered is used to plot and save violin plots.
#'
#' @param seur_obj Seurat object
#' @param marker_list gene name list that is to be plotted.
#' @param cols plot colour. Default is colour_palette() from the RNAseqEr library.
#' @param dir_lab a label that specifies which cell lineage is investigated. Needed
#' to structure output. THe default is "all_celltypes".
#' @param save_label Additional label that can be used to specify which condition
#' is used to save if a cluster or specific condition markers are tested. The default
#' is "cluster_mark" but it could be changed to "age", "treatment", "condition", etc.
#' which will be reflected in new directory names that save results.
#' @param split_by Secifies whether the data should be plotted seperately for a
#' specific metadata column.
#' @param save_dir Directory where the data should be saved.
#' @param numb_genes Number of genes that should be plotted together. Default is
#' 10.
#' @param plotheight height of the pdf output. Default is 20.
#' @param plotwidth width of the pdf output. Default is 8.
#' @param group_by Specifies which metadata column should be plotted on the x-axis
#' @param pt_size point size. Default is 0.1.
#' @param n_col number of columns in plot grit. Defaut is 1.
#'
#' @return creates a saved output to file.
#' @export
#'
#' @examples
#' markers <- c("SNAP25", "PLP1", "MAG", "CD74")
#' save_vln(cns, marker_list = markers)
#'
save_vln <- function(seur_obj,
                     marker_list,
                     cols = colour_palette(),
                     dir_lab = "all_celltypes",
                     save_label = "cluster_mark",
                     split_by = "NULL",
                     save_dir = getwd(),
                     numb_genes = 10,
                     plotheight = 20,
                     plotwidth = 8,
                     group_by = Idents(seur_obj),
                     pt_size = 0.1,
                     n_col = 1){
  gene_groups <- split(marker_list,
                       ceiling(seq_along(marker_list) / numb_genes))

  for(i in 1:length(gene_groups)){
    if(split_by != "NULL"){
      vln_plot <- VlnPlot(seur_obj,
                         features = unlist(gene_groups[i]),
                         group.by = group_by,
                         pt.size= pt_size,
                         ncol=1,
                         split.by = split_by,
                         cols = cols)
    }else{
      vln_plot <- VlnPlot(seur_obj,
                         features = unlist(gene_groups[i]),
                         group.by = group_by,
                         pt.size= pt_size,
                         ncol=1,
                         split.by = split_by,
                         cols = cols)
    }

    dir_save <- paste0(save_dir,
                       "/outs/",
                       dir_lab,
                       "/plots/",
                       save_label,
                       "/",
                       group_by,
                       "/violin_plots")
    if(dir.exists(dir_save) == FALSE){
      dir.create(dir_save, recursive = TRUE)
      print("creating new directory")
    }

    pdf(paste0(dir_save, "/", save_label, "_", i, ".pdf"),
        height=plotheight, width = plotwidth)
    print(vln_plot)
    dev.off()

    print(vln_plot)

  }
}

#' SCtype annotation
#' @param seur_obj Seurat object
#' @param custom_ref_genes
#' @param tissue_ref Specify which tissue is being used by chosing from the following:
#' "Brain", "Immune system", "Pancreas", "Liver", "Eye", "Kidney", "Brain", "Lung",
#' "Adrenal", "Heart", "Intestine", "Muscle", "Placenta", "Spleen", "Stomach",
#' "Thymus"
#' @param assay_use Which assay should be used for differential gene expression analysis.
#' Default is "RNA".
#' @param ident_use meta data column that contains clustering resolution of interest
#' for broad celltype annotation
#' @param low_confident_thres Relates celltype score to number of cells in that
#' category, the higher the value the more lenient. The default scTYpe threshold is 4.
#' If score is lower than the number of cells in that category divided by the
#' low_confident_thres, the cluster is annotated as "Unknown"
#' @param customclassif which name should be used for the meta data column saving the
#' scType annotation information. Default is "ScType_annotation".
#' @param save_dir Directory to be used for saving output.
#' @param dir_lab abel used for which data is analysed. Default is "all_celltypes"
#' @param plot_width width of DimPlot in output file. Default = 8.
#' @param plot_height height of DimPlot in output file. Default = 8.
#'
#' @return annotated seurat object
#' @import Seurat HGNChelper openxlsx
#' @export
#'
#' @examples
#' library(dplyr)
#' library(Seurat)
#' library(HGNChelper)
#' library(openxlsx)
#' data(cns)
#' cns <- sc_type_annot(cns,
#'                     ident_use = "rough_annot",
#'                     tissue_ref = "Brain",
#'                     save_plot = FALSE)
#'
#' print(names(cns@meta.data))
#'
sc_type_annot <- function(seur_obj,
                          custom_ref_genes = FALSE,
                          tissue_ref,
                          assay_use,
                          ident_use,
                          low_confident_thres = 4,
                          customclassif = "ScType_annotation",
                          save_plot = TRUE,
                          save_dir = getwd(),
                          dir_lab = "all_celltypes",
                          plot_width = 8,
                          plot_height = 8
                          ){

  ####
  #Load required data to run ScType, reconsider more elegant approach for loading
  # data

  lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = T)
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  if (custom_ref_genes == FALSE) {
    datbase <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypedatbasefull.xlsx"
    tissue_ref <- tissue_ref
  } else {
    datbase <- custom_gene_list
  }
  ####

  # Prepare ScType data

  gs_list <- gene_sets_prepare(datbase, tissue_ref)

  es_max <- sctype_score(
    scRNAseqData = seur_obj[[assay_use]]@scale.data, scaled = TRUE,
    gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
  )

  cL_resutls <- do.call(
    "rbind",
    lapply(
      unique(seur_obj@meta.data[[ident_use]]),
      function(cl) {
        es_max_cl <- sort(rowSums(es_max[,
                                         rownames(seur_obj@meta.data[seur_obj@meta.data[[ident_use]] ==
                                                                       cl, ])]),
                          decreasing = !0
        )
        head(data.frame(
          cluster = cl,
          type = names(es_max_cl),
          scores = es_max_cl,
          ncells = sum(seur_obj@meta.data[[ident_use]] == cl)
        ), 10)
      }
    )
  )


  sctype_scores <- cL_resutls %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)

  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) <
                       sctype_scores$ncells / low_confident_thres] <- "Unknown"

  #Annotate seurat object with scType annotation
  seur_obj@meta.data[[customclassif]] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    seur_obj@meta.data[[customclassif]][seur_obj@meta.data[[ident_use]] == j] <-
      as.character(cl_type$type[1])
  }

  if(save_plot == TRUE){
      dim_plot <- DimPlot(seur_obj,
                          reduction = plot_reduction,
                          label = plot_label,
                          repel = plot_repel,
                          group.by = customclassif
      )

      dim_plot_dir <- paste0(
        save_dir,
        "/outs/",
        dir_lab,
        "/plots/ScType_Annotated_plot"
      )

      if (dir.exists(dim_plot_dir) == FALSE) {
        dir.create(dim_plot_dir, recursive = TRUE)
        print("New directory created for saving DimPlot with ScType annotations.")
      }

      pdf(paste0(dim_plot_dir, "/ScType_annotated.pdf"),
          paper = "a4", width = plot_width, height = plot_height
      )

      print(dim_plot)

      dev.off()

      print(dim_plot)
  }

  return(seur_obj)

}


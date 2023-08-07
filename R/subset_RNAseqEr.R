#' subset_RNAseqEr
#' @description
#' This function subsets a Seurat object for each unique value in a provided
#' metadata column.
#'
#' @param seur_obj Seurat object
#' @param subset_column  column to use for subsetting
#' @param save_dir root directory in which outout is saved (subfolders are created)
#'
#' @return saves data to RDS files
#' @import here
#' @export
#'
#' @examples
#' library(here)
#' subset_RNAseqEr(cns,
#'                 subset_column = "rough_annot",
#'                 save_dir= here())
subset_RNAseqEr <- function(seur_obj,
                            subset_column = "RNAseqEr_annotation",
                            save_dir = getwd()){
  Idents(seur_obj) <- subset_column
  output_dir <- paste0(save_dir,
                       "/outs/data")
  dir.create(output_dir, recursive = TRUE)
  for(i in 1:length(levels(as.factor(seur_obj@meta.data[[subset_column]])))){
    curr_lev <- levels(as.factor(seur_obj@meta.data[[subset_column]]))[i]
    subs_seur <- subset(seur_obj, ident = curr_lev)
    save_name <- sub(" ", "_", curr_lev)
    saveRDS(subs_seur, paste0(output_dir, "/", save_name, ".RDS"))
  }
}

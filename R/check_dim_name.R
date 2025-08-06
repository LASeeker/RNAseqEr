#' Test name of dimesional reduction
#' @description
#' The function tests whether a provided name of a dimesional reduction is
#' present in its provided form in a Seurat object. It further test whether
#' the provided name is present in its letter sequence but converted
#' to upper or lower case . It then returns the provided plot reduction converted
#' to the version found in the seurat object which may be helpful particularly if
#' objects are converted between SingleCellExperiments and Seurat objects.
#' @param seur_obj Seurat object
#' @param plot_reduction name of dimensional reduction that is supposed to be
#' tested for its presence in the Seurat object and possibly adapted to upper/lower
#' case
#' @return plot_reduction if present in Seurat objet, ot its upper or lower case
#' version
#' @import dplyr
#' @export
#'
#' @examples
#' library(RNASeqEr)
#' library(dplyr)
#' data(cns)
#' dim_name <- check_dimred_name(seur_obj = cns,
#'                               plot_reduction = "umap")
#' print(dim_name)
#'
#' dim_name_2 <- check_dimred_name(seur_obj = cns,
#'                               plot_reduction = "UMAP")
#' print(dim_name_2)
#'
#' dim_name_3 <- check_dimred_name(seur_obj = cns,
#'                               plot_reduction = "TSNE")
#' print(dim_name_3)
#'
#' dim_name_4 <- check_dimred_name(seur_obj = cns,
#'                               plot_reduction = "NotPresent")
#' print(dim_name_4)
#'
check_dimred_name <- function(seur_obj,
                              plot_reduction = "umap"
                              ){
  red_test <- plot_reduction %in% names(seur_obj@reductions)
  if(red_test == FALSE){
    red_test2 = tolower(plot_reduction) %in% names(seur_obj@reductions)
    if(red_test2 == TRUE){
      plot_reduction = tolower(plot_reduction)
    }
    if(red_test2 == FALSE){
      red_test3 <- toupper(plot_reduction) %in% names(seur_obj@reductions)
      if(red_test3 == TRUE){
        plot_reduction = toupper(plot_reduction)
      }else{
        print("Check name of dimensional reduction. Not found in data")
        plot_reduction <- NA
      }
    }
  }
  return(plot_reduction)
}

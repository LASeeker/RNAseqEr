#' Create a list of resolutions that will be checked downstream for clustering
#' that can be validated
#'
#' @param seur_obj Seurat object
#' @param col_pattern Column pattern that should be used to recognise clustering
#' resolution. The default is "RNA_snn_res."
#'
#' @return a list of resolutions that will be used downstream
#' @export
#'
#' @examples
#' res_list <- select_res(cns, col_pattern = "RNA_snn_res.")
select_res <- function(seur_obj,
                       col_pattern = "RNA_snn_res."){
  extr_res_col <- grep(pattern = col_pattern, names(seur_obj@meta.data))
  res_names <- names(seur_obj@meta.data[extr_res_col])
  # gtools function, sorts gene_names alphanumeric:
  res_names <- gtools::mixedsort(res_names)
  for(i in 1: length(res_names)){
    curr_name <- res_names[i]
    curr_dat <- seur_obj@meta.data[[res_names[i]]]
    if(i == 1){
      res_df <- data.frame(matrix(1:ncol(seur_obj),
                                  nrow = ncol(seur_obj),
                                  ncol = 1))
      names(res_df) <- "row_number"
    }
    if(length(levels(as.factor(seur_obj@meta.data[[curr_name]]))) > 1){
      res_df[[curr_name]] <- curr_dat


    }
  }
  res_df <- res_df[, -1]

  #remove duplicated columns (resolution chage did not change clustering)
  res_df <- res_df[vapply(res_df, function(x) length(unique(x)) > 1, logical(1L))]

  # remove resolutions that have the same count of clusters as the previous
  level_df <- data.frame(matrix(NA, nrow= 1, ncol = 1))
  for(j in 1 : ncol(res_df)){
    level_df[,j] <- length(levels(as.factor((res_df[,j]))))
  }

  names(level_df) <- names(res_df)

  res_list <- list()
  res_list[1] <- names(level_df)[1]
  m <- 2

  for(k in 2: ncol(level_df)){
    if(level_df[,k - 1] != level_df[,k]){
      res_list[m] <- names(level_df)[k]
      m <- m + 1
    }

  }
  return(res_list)

}

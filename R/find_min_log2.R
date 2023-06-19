#' Find minimum log 2 fold change that is unlikely to be due to noise
#'
#' @description
#' Find a log2Fold threshold that for differential gene expression analysis which
#' is unlikely to be caused by noise in the dataset.
#' It is important to use a relatively homogeneous dataset (such as one cell type
#' of a control group only) as imput data. The function randomly subsets the
#' dataset of interest mixit times and performs differential gene expression anaysis.
#'
#' @param seur_obj Seurat object
#' @param maxit number of iterations run to find answer. Dataset is randomly
#' split into
#' @param test.to.use statistical test to use for differential gene expression
#' analysis. Default is "MAST", however, all other methods implemented in Seurat
#' work.
#' @param logfc_threshold mimimum log2fold threshold considered in differential
#' gene expression analysis to reduce run time. Default is set to a very low 0.05.
#' @param only_pos TRUE/FALSE whether positive or positive and negative log2fold
#' changes should be considered
#' @param divisor integral that specifies in how many parts the dataset should
#' be divided. Default is that the dataset is halved (2), however, if set to 3
#' for example two thrirds of the randomly subsetted data are compared with one
#' another.
#' @param percentile cutoff for percentile. Default is set to 0.95. If set to 1
#' it can be seen how large log2 fold changes can be if driven only by noise.
#'
#' @return returns a value which represents the percentile of log2 fold changes
#' that are driven by noise only.
#' @export
#'
#' @examples
#' library(Seurat)
#' min_log2 <- find_min_log2(pbmc_small, maxit = 10)
#' min_log2
#' summary(as.numeric(min_log2))
find_min_log2 <- function(seur_obj,
                          maxit = 100,
                          test.to.use = "MAST",
                          logfc_threshold = 0.05,
                          only_pos = TRUE,
                          divisor = 2,
                          percentile = 0.95) {
  my_list <- list()
  my_list_genes <- list()
  k <- 1

  for (i in 1:maxit) {
    curr_dat <- seur_obj
    max_count <- as.numeric(ncol(curr_dat))
    third_size <- ceiling(max_count / divisor)
    set.seed(k)
    sub_1 <- base::sample(
      x = colnames(curr_dat),
      size = third_size, replace = F
    )

    rest <- base::subset(
      x = colnames(curr_dat),
      subset = !(colnames(curr_dat) %in%
        sub_1)
    )

    sub_2 <- base::sample(
      x = rest,
      size = third_size - 1, replace = F
    )

    curr_dat$rep <- ifelse(rownames(curr_dat@meta.data) %in% sub_1,
      "rep_1",
      ifelse(rownames(curr_dat@meta.data) %in% sub_2,
        "rep_2", "NA"
      )
    )

    Idents(curr_dat) <- "rep"

    mark <- FindMarkers(curr_dat,
      ident.1 = "rep_1",
      ident.2 = "rep_2",
      test.use = test.to.use,
      logfc.threshold = logfc_threshold,
      only.pos = only_pos
    )
    if (nrow(mark) > 0) {
      perc <- quantile(mark$avg_log2FC, percentile)
    } else {
      perc <- 0
    }

    my_list[i] <- perc
    k <- k + 1
  }
  return(my_list)
}

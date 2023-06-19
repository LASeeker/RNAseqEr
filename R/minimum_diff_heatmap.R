#' Find meaningful differences in average gene expession between clusters
#' @description
#' Visual inspection of heatmaps may inform if clusters are too similar in their
#' gene expression and should be merged. When avoiding the visual evaluation and
#' comparing underpinning numeric differences, one has to establish first which
#' differences are meaningful. Here, a dataset that should be homogeneous (same
#' cell type and condition for example) is iteratively and randomly  split into
#' groups and the average gene expression across within both groups is compared.
#'
#' @param seur_obj seurat object that should be relatively uniform as in the same
#' celltyp from the same condition
#' @param maxit number of iterations used
#' @param divisor in how many groups is the dataset to be divided. Default is 2.
#' @param percentile Percentile used to find a cutoff for cluster similarity.
#' Default it 0.9 which means that across all iterations the maximum difference
#' between two randomly generated clusters and across all genes was in 90 % of the
#' iterations larger than value returned.
#' @param diff_threshold numeric threshold of similarity. Default is 10.
#'
#' @return A list of difference measures for each iteration
#' @export
#'
#' @examples
#' library(Seurat)
#' diff_val <- min_diff_hm(pbmc_small)
#' head(diff_val)
min_diff_hm <- function(seur_obj,
                        maxit = 100,
                        divisor = 2,
                        percentile = 0.9,
                        diff_threshold = 10){
  value_list <- list()
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

    curr_dat$grouping <- ifelse(rownames(curr_dat@meta.data) %in% sub_1,
                           "rep_1",
                           ifelse(rownames(curr_dat@meta.data) %in% sub_2,
                                  "rep_2", "NA"
                           )
    )



    cluster_averages_df <- as.data.frame(AverageExpression(curr_dat,
                                             group.by = "grouping",
                                             return.seurat = FALSE))
    colnames(cluster_averages_df) <- levels(as.factor(curr_dat@meta.data$grouping))


      max_diff <- abs(max(cluster_averages_df["rep_1"] - cluster_averages_df["rep_2"]))
      if(max_diff < diff_threshold){
        print(paste0("When considering iteration ",
                     i,
                     ", rep_1 and rep_2 should be merged as they appaear too similar with a value of ",
                     max_diff,
                     " ."))
      }else{
        print(max_diff)
      }

      value_list[i] <- max_diff

    k <- k + 1
  }
  val_ul <- unlist(value_list)
  print(summary(val_ul))
  perc <- quantile(val_ul, percentile)

  print(paste0("The ",
               percentile*100,
               "th percentile of differences between cluster equals ",
               round(perc, 2),
               "."))


  return(value_list)

  }


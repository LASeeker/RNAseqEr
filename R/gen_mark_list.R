#' gen)mark_list
#'
#' @param file_dir directory where differential gene expression results are saved
#' as .csv files.
#' @param ad_pval threshold for adjusted p-value, default is 0.05
#' @param avg_log threshold for minimum average log 2 fold change, default is 1.2
#' @param pct_1 threshold for percentage of cells/ nuclei within the cluster
#' of interest expressing the gene. Default is 0.25
#' @param pct_2 threshold for percentage of cells/ nuclei outside the cluster
#' of interest expressing the gene. Default is 0.6
#' @param pairwise TRUE/FALSE whether differential gene expression was performed
#' pairwise or not
#' @param n_top number of top genes per cluster considered. Default is 10.
#' @import Seurat dplyr here
#'
#' @return filtered data frame of potential cluster marker genes with n_top candidate
#' genes per cluster/ group of interest
#' @export
#'
#' @examples
#' library(Seurat)
#' library(dplyr)
#' library(here)
#' save_dir_cr <- here("data", "dge")
#' dir.create(save_dir_cr, recursive = TRUE)
#' df_0_8 <- as.data.frame(dge_0_8)
#' df_1_0 <- as.data.frame(dge_1_0)
#' write.csv(df_0_8, here(save_dir_cr, "dge_0_8.csv"))
#' write.csv(df_1_0, here(save_dir_cr, "dge_1_0.csv"))
#' clu_mark <- gen_mark_list(file_dir = here("data", "dge"))
gen_mark_list <- function(file_dir = getwd(),
                          ad_pval = 0.05,
                          avg_log = 1.2,
                          pct_1 = 0.25,
                          pct_2 = 0.6,
                          pairwise = FALSE,
                          n_top = 10) {
  temp <- list.files(
    path = file_dir,
    pattern = "*.csv"
  )
  myfiles <- lapply(paste(file_dir, temp, sep = "/"), read.csv)

  for (i in 1:length(myfiles)) {
    dat <- myfiles[[i]]
    av_log_fil <- subset(dat, dat$avg_log2FC > avg_log &
        dat$pct.1 > pct_1 &
        dat$pct.2 < pct_2 &
        dat$p_val_adj < ad_pval)
    if (pairwise == TRUE) {
        av_log_fil <- av_log_fil[order("avg_log2FC"), ]
        top10 <- av_log_fil %>% top_n(n_top, avg_log2FC)
        tail10 <- tail(av_log_fil, n_top)
        top10$gene <- top10$X
        tail10$gene <- tail10$X
        top10_df <- rbind(top10, tail10)
    } else {
        av_log_fil$cluster <- as.character(av_log_fil$cluster)
        top10 <- av_log_fil %>%
          group_by(cluster) %>%
          top_n(n_top, avg_log2FC)
      top10_df <- as.data.frame(top10)
      }

      if (exists("fil_genes") == FALSE) {
        fil_genes <- top10_df
      } else {
        fil_genes <- rbind(fil_genes, top10_df)
      }


    fil_genes <- fil_genes[!duplicated(fil_genes$gene), ]
  }

  return(fil_genes)
}

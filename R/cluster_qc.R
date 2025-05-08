#' cluster_qc
#' @description
#' This function performs a quick cluster quality control where it checks it certain
#' levels of metadata factors contribute more to certain clusters than others.
#' It can be used for example to check if certain donor_ids or sample_ids
#' drive clustering.
#' @param seur_obj Seurat objects which is pre-processed and clustered
#' @param vars this is a string of metadata column names that is to be checked.
#' it can be "sample_id" or c("sample_id", "donor_id" for example.)
#' @param cluster_col column in metadata that contains cluster resolution that is
#' to be checked
#' @param thres threshold of which proportion of vars should ideally contribute
#' to each cluster. The default is set to 0.7 but may need to be adjusted based
#' on the experimental design.
#' @param low_thres lower threshold of proportion of vars that need to be exceeded
#' for cluster to pass QC. Default is 0.1.
#' @param colors plot colors. Default os colour_palette() from RNAseqEr library.
#' @param save_dir diretory in which output should be saved
#' @param dir_lab label of current cell type studies, default is "all_celltypes".
#' @param width width of output plots. Default is 8.
#' @param height height of output plots tefault os 5.
#' @param scale  scale used to flag clusters. The smaller the more stringend is
#' cluster QC. Default is 15.
#'
#' @return seurat object with updated cluster QC measures. Additionally, plots
#' and tables are saved to file
#' @import Seurat
#' @export
#'
#' @examples
#' seur_obj <- cluster_qc(cns,
#'                    cluster_col = "rough_annot",
#'                    vars = c("caseNO", "process_number"))
cluster_qc <- function(seur_obj,
                       vars,
                       cluster_col,
                       thres = 0.7,
                       low_thres = 0.1,
                       colors = colour_palette(),
                       save_dir= getwd(),
                       dir_lab = "all_celltypes",
                       width=8,
                       height= 5,
                       scale = 15
                       ){
  # count how many cells there are in each group and cluster
  met_dat <- seur_obj@meta.data

  for(i in 1: length(vars)){
      curr_factor <- as.factor(met_dat[[vars[i]]])

      total <- length(levels(curr_factor))

      threshold <- total * thres
      low_threshold <- total * low_thres


      low_threshold <- ifelse(low_threshold > 2, low_threshold, 2)


      curr_tab <- table(met_dat[[cluster_col]] , met_dat[[vars[i]]])
      # For each cluster (on the rows) sum of individuals that do have cells on that cluster
      curr_df <- as.data.frame(rowSums(curr_tab > 0))
      names(curr_df) <- paste0(vars[i], "_per_cluster")

      save_df <- data.frame(cluster = rownames(curr_df),
                            X = curr_df[[paste0(vars[i], "_per_cluster")]])

      save_df$threshold<- threshold
      save_df$threshold_low <- low_threshold

      save_df$placeholder <- ifelse(save_df$X < threshold,
                                    paste0(" Less ",
                                           vars[i],
                                           " contribute to cluster than set threshold"
                                           ),
                                    paste0("More than ",
                                           thres*100,
                                           " % of ",
                                           vars[i],
                                           " contribute to this cluster."))

      save_df$placeholder_2 <- ifelse(save_df$X < low_threshold,
                                    paste0("* Less than than ", low_thres *100,
                                           " % of ",
                                           vars[i],
                                           " contribute to this cluster."),
                                    paste0("* QC pass"))
      save_df$placeholder_3 <- ifelse(save_df$X < low_thres,
                                      paste0("Remove"),
                                      ifelse(save_df$X > low_thres &
                                               save_df$X < threshold,
                                             paste0("Watch ", vars[i], " variation"),
                                             paste0("Pass")))
      names(save_df)[1] <- paste0(cluster_col)

      names(save_df)[ncol(save_df)-5] <- paste0(vars[i], "_per_cluster")


      names(save_df)[ncol(save_df)-2] <- paste0(vars[i], "_per_cluster_qc")

      names(save_df)[ncol(save_df)-1] <- paste0(vars[i], "_per_cluster_qc_low")

      names(save_df)[ncol(save_df)] <- paste0(vars[i], "_conclusion")



      # Check for individual entries contributing more than a set proportion to
      # clusters

      prop_var_table <- prop.table(curr_tab, margin = 1)
      # change the format to be a data.frame, this also expands to long formatting
      prop_var_df <- as.data.frame(prop_var_table)
      colnames(prop_var_df) <- c("cluster", paste0(vars[i]), "proportion")
      flag_clu <- prop_var_df[prop_var_df$proportion > (100/total)/scale, ]

      save_df$Y <- ifelse(save_df[[cluster_col]] %in% flag_clu[[cluster_col]],
                          paste0("* Individual entities of ",
                                 vars[i],
                                 " may be proportionally overrepresented in this cluster."),
                          paste0(" * No individual ",
                                 vars[i],
                                 " contributed more than ",
                                 threshold[i],
                                 " to this cluster. (QC pass)"))


      names(save_df)[ncol(save_df)] <- paste0(vars[i], "_proportional_contrib")


      curr_plot <- ggplot(data = prop_var_df,
                           aes(x = cluster, y = proportion, fill = .data[[vars[i]]] )) +
        geom_bar(stat = "identity") + theme_classic() +
        scale_fill_manual(values = colors) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      curr_plot

      plot_dir <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/plots/cluster_qc"
                         )

      if(dir.exists(plot_dir) == FALSE){
          dir.create(plot_dir)
      }


      pdf(paste0(plot_dir,
                 "/",
                 vars[i], "_cluster_qc.pdf"), width=width, height=height)
      plot(curr_plot)

      dev.off()

      plot(curr_plot)


      tab_dir <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/tables/cluster_qc"
      )
      if(dir.exists(tab_dir) == FALSE){
        dir.create(tab_dir)
      }

      write.csv(save_df, paste0(tab_dir, "/", vars[i], "_cluster_qc_tab.csv"))


      # Add summary QC info to seurat object
      new_dat <- merge(met_dat, save_df, by = cluster_col, all = TRUE)

      order <- match(met_dat[[cluster_col]], new_dat[[cluster_col]])

      new_dat <- new_dat[order, ]

      new_col_name <- paste0("RNAseqEr_clust_qc_", vars[i])

      seur_obj@meta.data[[new_col_name]] <- new_dat[[ncol(new_dat)-1]]

      dim_plot <- DimPlot(seur_obj,
              group.by = new_col_name,
              cols = colour_palette())


      plot_dir_2 <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/plots/cluster_qc/DimPlots"
      )

      if(dir.exists(plot_dir_2) == FALSE){
        dir.create(plot_dir_2)
      }


      pdf(paste0(plot_dir_2,
                 "/",
                 vars[i], "_cluster_qc_dimPlot.pdf"), width=width, height=height)
      plot(dim_plot)

      dev.off()

      plot(dim_plot)


  }
  return(seur_obj)


}














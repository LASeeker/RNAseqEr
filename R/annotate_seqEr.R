#' Automatic annotation of broad celltypes using ScType
#' @description
#' ScType (https://www.nature.com/articles/s41467-022-28803-w) is used to identify
#' cell types automatically. If a onsiderable number of cells remain unknown,
#' an additional step which is based on differentially expressed genes is used
#' to attempt an celltype identification.
#'
#' @param seur_obj Seurat object
#' @param tissue_ref Specify which tissue is being used by chosing from the following:
#' "Brain", "Immune system", "Pancreas", "Liver", "Eye", "Kidney", "Brain", "Lung",
#' "Adrenal", "Heart", "Intestine", "Muscle", "Placenta", "Spleen", "Stomach",
#' "Thymus"
#' @param ident_use meta data column that contains clustering resolution of interest
#' for broad celltype annotation
#' @param test_use statistical test used for differential gene expression analysis.
#' Default is "MAST".
#' @param min_pct minimum percentage of cells within the cluster of interest
#' expressing the gene. Default is 0.25.
#' @param logfc_threshold Minimal log2 fold change for differential gene expression.
#' Default is 0.8.
#' @param assay_use Which assay should be used for differential gene expression analusis.
#' Default is "RNA".
#' @param save_dir Directory to be used for saving output.
#' @param custom_ref_genes TRUE/FALSE whether a list of genes and celltypes is
#' being provided by the user.
#' @param custom_gene_list if custom_ref_genes == TRUE provide custom gene list here.
#' @param customclassif which name should be used for the meta data column saving the
#' scType annotation information. Default is ScType_annotation"
#' @param plot_reduction which dimensional reduction should be used to plot data
#' with ScType labels. Default is "umap"
#' @param plot_label TRUE/FALSE whether ScType/ RNAseqEr annotation  labels should
#' be displayed on DimPlots. Default is TRUE.
#' @param plot_repel TRUE/FALSE whether ScType/ RNAseqEr annotation labels should
#' be repelled on DimPlots. Default is TRUE.
#' @param plot_width width of DimPlot in output file. Default = 8.
#' @param plot_height height of DimPlot in output file. Default = 8.
#' @param dge_present TRUE/FALSE if previously generated DGE results are present.
#' Should be set to FALSE at each new tested ident_use for the first time, but
#' can be set to TRUE thereafter to speed up process.
#' @param dge TURE/FALSE whether differential gene expression should be performed.
#' Default is TRUE, to allow for annotation of celltypes that were not annotated
#' based on scType alone.
#' @param proportion threshold proportion for annotation. Default is set to 3 which
#' means that the number differentially expressed genes in a cluster must be
#' larger than a third of the total number of genes that define a cell type in scTyoe to
#' be cosidered that cell type. Increasing the number makes annotation less stringent,
#' decreasing it makes it more stringent.
#' @param colours Colours for DimPlots. Default is colour_pelette() from The RNAseqEr
#' library
#'
#' @return Seurat object with added ScType and RNAseqEr annotation
#' @export
#' @import dplyr Seurat HGNChelper openxlsx
#'
#' @examples
#' library(dplyr)
#' library(Seurat)
#' library(HGNChelper)
#' library(openxlsx)
#' data(cns)
#' cns <- annotate_seqEr(cns, ident_use = "rough_annot")
#'
#' # or:
#' data(cns)
#' cns <- annotate_seqEr(cns, ident_use = "Fine_cluster")
annotate_seqEr <- function(seur_obj,
                           ident_use,
                           tissue_ref = "Brain",
                           test_use = "MAST",
                           min_pct = 0.25,
                           logfc_threshold = 0.8,
                           assay_use = "RNA",
                           save_dir = getwd(),
                           custom_ref_genes = FALSE,
                           custom_gene_list,
                           customclassif = "ScType_annotation",
                           plot_reduction = "umap",
                           plot_label = TRUE,
                           plot_repel = TRUE,
                           plot_width = 8,
                           plot_height = 8,
                           dge_present = FALSE,
                           dge = TRUE,
                           proportion = 3,
                           colours = colour_palette()) {
  lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = T)
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  if (custom_ref_genes == FALSE) {
    db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    tissue_ref <- tissue_ref
  } else {
    db_ <- custom_gene_list
  }


  gs_list <- gene_sets_prepare(db_, tissue_ref)

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
    sctype_scores$ncells / 4] <- "Unknown"
  print(sctype_scores[, 1:3])


  seur_obj@meta.data[[customclassif]] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    seur_obj@meta.data[[customclassif]][seur_obj@meta.data[[ident_use]] == j] <-
      as.character(cl_type$type[1])
  }


  dim_plot <- DimPlot(seur_obj,
    reduction = plot_reduction,
    label = plot_label,
    repel = plot_repel,
    group.by = customclassif
  )

  dim_plot_dir <- paste0(
    save_dir,
    "/outs/plots/ScType_Annotated_plot"
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

  marker_dir <- paste0(
    save_dir,
    "/outs/tables/broad_celltype_markers"
  )

  if (dir.exists(marker_dir) == FALSE) {
    dir.create(marker_dir, recursive = TRUE)
    print("New directory created for saving differential gene expression
                results for broad clusters that are assumed to represent celltypes.")
  }



  if(dge_present == TRUE){

    markers <- read.csv(paste0(marker_dir,
                               "/ScType_cluster_markers_",
                               ident_use,
                               ".csv"))
    #remove previous annotation if present
    bool <- grepl("new_annot",
                  names(markers), fixed = FALSE)
    markers <- markers[,!bool]
    bool <- grepl("RNAseqEr_ann",
                  names(markers), fixed = FALSE)
    markers <- markers[,!bool]
    bool <- grepl("X.",
                  names(markers), fixed = FALSE)
    markers <- markers[,!bool]

  }


  Idents(seur_obj) <- customclassif
  seur_obj$scType_prelim_ann <- ifelse(seur_obj@meta.data[[customclassif]] == "Unknown",
                                       paste(seur_obj@meta.data[[ident_use]],
                                             seur_obj@meta.data[[customclassif]],
                                             sep = "_"
                                       ),
                                       paste(seur_obj@meta.data[[customclassif]])
  )
  Idents(seur_obj) <- "scType_prelim_ann"



  if (dge == TRUE) {

    markers <- FindAllMarkers(seur_obj,
        min.pct = min_pct,
        test.use = test_use,
        assay = assay_use,
        only.pos = TRUE,
        logfc.threshold = logfc_threshold
      )
  }




  levels <- levels(as.factor(markers$cluster))




  if ("Unknown" %in% levels(as.factor(seur_obj@meta.data[[customclassif]]))) {

    level_list <- list()
      for (g in 1:length(levels)) {
        level_list[g] <- grepl("Unknown", levels[g], fixed = FALSE)
      }


      unknown_levels <- levels[unlist(level_list)]


      if (length(unknown_levels > 0)) {
        for (o in 1:length(unknown_levels)) {
          temp_dat <- subset(markers, markers$cluster == unknown_levels[o])
          fil_temp <- subset(
            temp_dat,
            temp_dat$p_val_adj < 0.05 &
              temp_dat$avg_log2FC > 1.5
          )

          cell_type_list <- names(gs_list$gs_positive)
          for (l in 1:length(cell_type_list)) {
            curr_ct <- cell_type_list[l]
            bool <- fil_temp$gene %in% gs_list$gs_positive[[curr_ct]]
            fil_mark <- fil_temp[bool, ]
            if (nrow(fil_mark) > length(gs_list$gs_positive[[curr_ct]]) / proportion) {
              if (exists("df_annot") == FALSE) {
                df_annot <- data.frame(
                  cluster = unknown_levels[o],
                  new_annot = curr_ct
                )
              } else {
                df_prelim <- data.frame(
                  cluster = unknown_levels[o],
                  new_annot = curr_ct
                )
                df_annot <- rbind(df_annot, df_prelim)
                print(l)
              }
            }
          }
          markers_add <- merge(markers, df_annot,
            by = "cluster",
            all = TRUE
          )


          markers_add$RNAseqEr_ann <- ifelse(markers_add$cluster %in% unknown_levels,
            paste(markers_add$new_annot),
            paste(markers_add$cluster)
          )


          write.csv(markers_add, paste0(marker_dir,
                                        "/ScType_cluster_markers_",
                                        ident_use,
                                        ".csv"))


          df_save <- data.frame(
            scType_prelim_ann = markers_add$cluster,
            RNAseqEr_annotation = markers_add$RNAseqEr_ann
          )
          df_save_uniq <- unique(df_save)

          met_dat <- seur_obj@meta.data



          new_met_dat <- merge(met_dat,
            df_save_uniq,
            by ="scType_prelim_ann",
            all = TRUE
          )


          print(head(new_met_dat))

          new_met_dat <- new_met_dat[match(
            met_dat$scType_prelim_ann,
            new_met_dat$scType_prelim_ann
          ), ]


          bool <- grepl("RNAseqEr_annotation", names(new_met_dat))
          ann_use <- as.data.frame(new_met_dat[, bool])


         if(ncol(ann_use) > 1){
            ann_use <- as.data.frame(ann_use[,2])
         }
          names(ann_use) <- "RNAseqEr_annotation"


          #seur_obj@meta.data$RNAseqEr_annotation <- new_met_dat$RNAseqEr_annotation
          seur_obj@meta.data$RNAseqEr_annotation <- ann_use$RNAseqEr_annotation

        }
        }
      }else{
        print("debug horse")
          seur_obj@meta.data$RNAseqEr_annotation <- seur_obj@meta.data[[customclassif]]
          write.csv(markers, paste0(marker_dir,
                                        "/ScType_cluster_markers_",
                                        ident_use,
                                        ".csv"))
      }


        dim_plot <- DimPlot(seur_obj,
                              group.by = "RNAseqEr_annotation",
                              cols = colours,
                              label = TRUE)
         print(dim_plot)

         pdf(paste0(dim_plot_dir, "/RNAseqEr_annotated.pdf"),
             paper = "a4", width = plot_width, height = plot_height
         )

         print(dim_plot)

         dev.off()


         if("Unknown" %in% levels(as.factor(seur_obj@meta.data[[customclassif]]))){
           if(exists("unknown_levels") == TRUE){
             if(length(unknown_levels) > 1){
              seur_obj@meta.data$RNAseqEr_annotation_detailed <-
                ifelse(seur_obj@meta.data$ScType_annotation == "Unknown",
                    paste(seur_obj@meta.data$RNAseqEr_annotation,
                          seur_obj@meta.data[[ident_use]],
                          sep = "_"),
                    paste(seur_obj@meta.data$RNAseqEr_annotation))

              dim_plot <- DimPlot(seur_obj,
                               group.by = "RNAseqEr_annotation_detailed",
                               cols = colours,
                               label = TRUE)
           print(dim_plot)

           pdf(paste0(dim_plot_dir, "/RNAseqEr_annotated_detailed.pdf"),
               paper = "a4", width = plot_width, height = plot_height)
           print(dim_plot)

           dev.off()
           }
           }
         }
         return(seur_obj)
      }



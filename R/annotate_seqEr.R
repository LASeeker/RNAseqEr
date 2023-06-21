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
#' @param ident meta data column that contains clustering resolution of interest
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
#' @param cusom_gene_list if custom_ref_genes == TRUE provide custom gene list here.
#' @param customclassif which name should be used for the meta data column saving the
#' scType annotation information. Default is ScType_annotation"
#' @param plot_reduction which dimensional reduction should be used to plot data
#' with ScType labels. Default is "umap"
#' @param plot_label TRUE/FALSE whether ScType labels should be displayed on the
#' DimPlot. Default is TRUE.
#' @param plot_repel TRUE/FALSE whether ScType labels should be repelled on the
#' DimPlot. Default is TRUE.
#' @param plot_width width of DimPlot in output file. Default = 8.
#' @param plot_height height of DimPlot in output file. Default = 8.
#' @param dge TURE/FALSE whether differential gene expression should be perfomed.
#' Default is FALSE, because it may take some time.
#'
#' @return seurat object with added ScType annotation
#' @export
#' @import dplyr Seurat HGNChelper openxlsx
#'
#' @examples
#' library(dplyr)
#' library(Seurat)
#' library(HGNChelper)
#' library(openxlsx)
#' pbmc_small <- annotate_seqEr(pbmc_small,
#'     tissue_ref= "Immune system",
#'     ident = "RNA_snn_res.1")
annotate_seqEr <- function(seur_obj,
                           ident,
                           tissue_ref = "Brain",
                           test_use = "MAST",
                           min_pct = 0.25,
                           logfc_threshold = 0.8,
                           assay_use = "RNA",
                           save_dir = get_wd(),
                           custom_ref_genes = FALSE,
                           cusom_gene_list,
                           customclassif = "ScType_annotation",
                           plot_reduction = "umap",
                           plot_label = TRUE,
                           plot_repel = TRUE,
                           plot_width = 8,
                           plot_height = 8,
                           dge = FALSE){
    lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R");
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

    # load gene set preparation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    # load cell type annotation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

    if(custom_ref_genes == FALSE){
    db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";tissue = tissue_ref
    }else{
      db_ <- cusom_gene_list
    }

    gs_list = gene_sets_prepare(db_, tissue)

    es_max = sctype_score(scRNAseqData = seur_obj[[assay]]@scale.data, scaled = TRUE,
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


    cL_resutls = do.call("rbind",
                         lapply(unique(seur_obj@meta.data[[ident]]),
                                         function(cl){
      es_max_cl = sort(rowSums(es_max[ ,rownames(seur_obj@meta.data[seur_obj@meta.data[[ident]]==cl, ])]),
                       decreasing = !0)
      head(data.frame(cluster = cl,
                      type = names(es_max_cl),
                      scores = es_max_cl,
                      ncells = sum(seur_obj@meta.data[[ident]]==cl)), 10)
    }))

    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) <
                         sctype_scores$ncells/4] = "Unknown"
    print(sctype_scores[,1:3])


    seur_obj@meta.data[[customclassif]] = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,];
      seur_obj@meta.data[[customclassif]][seur_obj@meta.data[[ident]] == j] =
        as.character(cl_type$type[1])
    }

    dim_plot <- DimPlot(seur_obj,
                        reduction = plot_reduction,
                        label = plot_label,
                        repel = plot_repel,
                        group.by = customclassif)

    dim_plot_dir <- paste0(save_dir,
                         "/outs/plots/ScType_Annotated_plot"
    )

    if(dir.exists(dim_plot_dir) == FALSE){
      dir.create(dim_plot_dir, recursive = TRUE)
      print("New directory created for saving DimPlot with ScType annotations.")
    }

    pdf(paste0(dim_plot_dir, "/ScType_annotated.pdf"),
        paper="a4", width=plot_width, height = plot_height)

    print(dim_plot)

    dev.off()

    print(dim_plot)


    if(dge == TRUE){
        Idents(seur_obj) <- customclassif


        markers <- FindAllMarkers(seur_obj,
                                  min.pct = min_pct,
                                  test.use = test_use,
                                  assay = assay_use,
                                  only.pos = TRUE,
                                  logfc.threshold = logfc_threshold)

        marker_dir <- paste0(save_dir,
                             "/outs/tables/broad_celltype_merkers"
                             )

        if(dir.exists(marker_dir) == FALSE){
          dir.create(marker_dir, recursive = TRUE)
          print("New directory created for saving differential gene expression
                results for broad clusters that are assumed to represent celltypes.")
        }


        levels <- levels(as.factor(markers$cluster))

        if("Unknown" %in% levels){
            temp_dat <- subset(markers, markers$cluster == "Unknown")
            fil_temp <- subset(temp_dat,
                               temp_dat$p_val_adj < 0.05 &
                                 temp_dat$avg_log2FC > 1.5)

            cell_type_list <- names(gs_list$gs_positive)
            for(l in 1: length(cell_type_list)){
              curr_ct <- cell_type_list[l]
              bool <- fil_temp$gene %in% gs_list$gs_positive[[curr_ct]]

              fil_mark <- fil_temp[bool,]
              if(nrow(fil_mark) > length(gs_list$gs_positive[[curr_ct]])/3){
                markers$annotated <- ifelse(markers$cluster == "Unknown",
                                             curr_ct,
                                             paste(markers$cluster))
                seur_obj@meta.data$RNAseqEr_annot <-
                  ifelse(seur_obj@meta.data[[customclassif]] == "Unknown",
                         paste(markers$annotated)[1],
                         paste(seur_obj@meta.data[[customclassif]]))
                }

          }
        }

        write.csv(markers, paste0(marker_dir, "/ScType_cluster_markers.csv"))
    }

    return(seur_obj)

    }


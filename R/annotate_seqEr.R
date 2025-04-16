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
#' @param dir_lab label used for which data is analysed. Default is "all_celltypes"
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
                           dir_lab = "all_celltypes",
                           custom_ref_genes = TRUE,
                           custom_gene_list = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",
                           customclassif = "ScType_annotation",
                           plot_reduction = "umap",
                           plot_label = TRUE,
                           plot_repel = TRUE,
                           plot_width = 8,
                           plot_height = 8,
                           dge = TRUE,
                           proportion = 3,
                           unknown_threshold = 0.25, # New parameter for controlling what's considered "Unknown"
                           colours = colour_palette()) {
  
  # Load required ScType functions
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  # Create directories for output
  dim_plot_dir <- paste0(save_dir, "/outs/", dir_lab, "/plots/ScType_Annotated_plot")
  marker_dir <- paste0(save_dir, "/outs/", dir_lab, "/tables/DGE/broad_celltype_markers")
  
  if (!dir.exists(dim_plot_dir)) {
    dir.create(dim_plot_dir, recursive = TRUE)
    print("New directory created for saving DimPlot with ScType annotations.")
  }
  
  if (!dir.exists(marker_dir)) {
    dir.create(marker_dir, recursive = TRUE)
    print("New directory created for saving differential gene expression results.")
  }
  
  # Step 1: Prepare gene sets
  if (custom_ref_genes == FALSE) {
    db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  } else {
    db_ <- custom_gene_list
  }
  
  gs_list <- gene_sets_prepare(db_, tissue_ref)
  
  # Step 2: Initial ScType annotation on expressed genes
  es_max <- sctype_score(
    scRNAseqData = seur_obj[[assay_use]]@scale.data, 
    scaled = TRUE,
    gs = gs_list$gs_positive, 
    gs2 = gs_list$gs_negative
  )
  
  # Assign cell types based on clusters
  cL_results <- do.call(
    "rbind",
    lapply(
      unique(seur_obj@meta.data[[ident_use]]),
      function(cl) {
        cells_in_cluster <- rownames(seur_obj@meta.data[seur_obj@meta.data[[ident_use]] == cl, ])
        es_max_cl <- sort(rowSums(es_max[, cells_in_cluster]), decreasing = TRUE)
        head(data.frame(
          cluster = cl,
          type = names(es_max_cl),
          scores = es_max_cl,
          ncells = length(cells_in_cluster)
        ), 10)
      }
    )
  )
  
  # Get top score for each cluster
  sctype_scores <- cL_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)
  
  # Mark clusters with low confidence scores as "Unknown"
  # Now using the user-configurable unknown_threshold parameter
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < 
                       sctype_scores$ncells * unknown_threshold] <- "Unknown"
  print(sctype_scores[, 1:3])
  
  # Add ScType annotation to Seurat object
  seur_obj@meta.data[[customclassif]] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    seur_obj@meta.data[[customclassif]][seur_obj@meta.data[[ident_use]] == j] <-
      as.character(cl_type$type[1])
  }
  
  # Save initial ScType plot
  p1 <- DimPlot(seur_obj,
                reduction = plot_reduction,
                label = plot_label,
                repel = plot_repel,
                group.by = customclassif)
  
  pdf(paste0(dim_plot_dir, "/ScType_annotated.pdf"),
      paper = "a4", width = plot_width, height = plot_height)
  print(p1)
  dev.off()
  
  # Step 3: Perform DE analysis if requested
  if (dge) {
    # Set identity for finding markers
    Idents(seur_obj) <- ident_use
    
    # Find markers for all clusters
    all_markers <- FindAllMarkers(seur_obj,
                                  min.pct = min_pct,
                                  test.use = test_use,
                                  assay = assay_use,
                                  only.pos = TRUE,
                                  logfc.threshold = logfc_threshold)
    
    # Save all markers
    write.csv(all_markers, paste0(marker_dir, "/All_cluster_markers_", ident_use, ".csv"))
    
    # Create preliminary annotation column combining cluster ID and ScType annotation
    seur_obj$scType_prelim_ann <- ifelse(seur_obj@meta.data[[customclassif]] == "Unknown",
                                         paste(seur_obj@meta.data[[ident_use]],
                                               seur_obj@meta.data[[customclassif]],
                                               sep = "_"),
                                         paste(seur_obj@meta.data[[customclassif]]))
    
    # Step 4: Enhanced annotation of "Unknown" clusters using DE genes
    unknown_clusters <- unique(seur_obj@meta.data[[ident_use]][seur_obj@meta.data[[customclassif]] == "Unknown"])
    
    if (length(unknown_clusters) > 0) {
      # Initialize final annotation column
      seur_obj$RNAseqEr_annotation <- seur_obj@meta.data[[customclassif]]
      
      # Create dataframe for storing new annotations
      annotation_df <- data.frame(cluster = character(), new_annotation = character())
      
      # Process each unknown cluster
      for (cluster_id in unknown_clusters) {
        # Get DE genes for this cluster
        cluster_markers <- subset(all_markers, all_markers$cluster == cluster_id & 
                                    all_markers$p_val_adj < 0.05 & 
                                    all_markers$avg_log2FC > 1.5)
        
        best_match <- NULL
        best_match_count <- 0
        
        # Check against each cell type in our reference
        for (cell_type in names(gs_list$gs_positive)) {
          # How many marker genes match this cell type's reference genes
          matches <- sum(cluster_markers$gene %in% gs_list$gs_positive[[cell_type]])
          required_matches <- length(gs_list$gs_positive[[cell_type]]) / proportion
          
          # If we have enough matches and it's better than previous best
          if (matches > required_matches && matches > best_match_count) {
            best_match <- cell_type
            best_match_count <- matches
          }
        }
        
        # If we found a match, record it
        if (!is.null(best_match)) {
          annotation_df <- rbind(annotation_df, 
                                 data.frame(cluster = cluster_id, 
                                            new_annotation = best_match))
        }
      }
      
      # Apply new annotations to the Seurat object
      if (nrow(annotation_df) > 0) {
        for (i in 1:nrow(annotation_df)) {
          cluster_id <- annotation_df$cluster[i]
          new_label <- annotation_df$new_annotation[i]
          
          # Update the annotation for this cluster
          seur_obj$RNAseqEr_annotation[seur_obj@meta.data[[ident_use]] == cluster_id] <- new_label
        }
        
        # Save annotation details
        write.csv(annotation_df, paste0(marker_dir, "/Unknown_cluster_new_annotations_", ident_use, ".csv"))
      }
    }
    
    # Create detailed annotation that includes cluster ID for remaining unknowns
    if ("Unknown" %in% levels(as.factor(seur_obj$RNAseqEr_annotation))) {
      seur_obj$RNAseqEr_annotation_detailed <- ifelse(
        seur_obj$RNAseqEr_annotation == "Unknown",
        paste(seur_obj$RNAseqEr_annotation, seur_obj@meta.data[[ident_use]], sep = "_"),
        paste(seur_obj$RNAseqEr_annotation)
      )
      
      # Plot detailed annotations
      p3 <- DimPlot(seur_obj,
                    group.by = "RNAseqEr_annotation_detailed",
                    cols = colours,
                    label = TRUE)
      
      pdf(paste0(dim_plot_dir, "/RNAseqEr_annotated_detailed.pdf"),
          paper = "a4", width = plot_width, height = plot_height)
      print(p3)
      dev.off()
    }
    
    # Plot final annotations
    p2 <- DimPlot(seur_obj,
                  group.by = "RNAseqEr_annotation",
                  cols = colours,
                  label = TRUE)
    
    pdf(paste0(dim_plot_dir, "/RNAseqEr_annotated.pdf"),
        paper = "a4", width = plot_width, height = plot_height)
    print(p2)
    dev.off()
  }
  
  return(seur_obj)
}



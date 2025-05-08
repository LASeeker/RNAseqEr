#' Fast Seurat Processing of RNAseq data until dimensional reduction
#' @description This function performs the standard Seurat processing on a quality
#' controlled Seurat object an returns a Seurat object with added dimensional reduction
#' and clustering at multiple resolutions. Now works for sketched data, too.
#'
#' @param seur_obj A Seurat object that has been quality controlled and if
#' necessary batch corrected and integrated
#' @param n_pcs number of dimensions (1:n_pcs) used for findNeighbors() and and
#' non linear dimensional reduction (UMAP, TSNE).
#'
#' @param res resolutions for clustering
#' @param select_genes genes that should be scaled. Default is all genes in Seurat
#' object but setting to variable genes may be valid too.
#' @param elbow_dims How many dimensions should be plotted in elbow plot to
#' determine n_pcs
#' @param tsne TRUE/FALSE, default is true which means TSNE dim reduction is
#' performed additionally to UMAP
#' @param dir_lab label that indicated cell lineage to structure output. The
#' default is "all_celltypes"
#' @param save_dir directory where the output should be saved. Default is the
#' current working directory.
#' @param plotheight height of elbow plot output. Default is 5.
#' @param plotwidth width of elbow plot output. Default is 5.
#' @param sketch whether data should be sketched (see Seurat documentation), useful
#' for very large datasets. Default is FALSE.
#' @param sketch_ncells number of cells used for sketched dataset. Default is 50000
#' @param sketch_method  method used for sketching. Default is "LeverageScore", 
#' alternative is "Uniform".
#' @param sketched_assay_name name of sketched assay output, Default is "sketch".
#' @return returns a Seeurat object that is normalised, has information on linear and
#' non-linear dimensional reductions (PCA, UMAP, TSNE) and is clustered at different
#' resolutions.
#' @import Seurat
#' @export
#'
#' @examples
#' library(Seurat)
#' cns <- seurat_proc(cns)
#'
seurat_proc <- function(seur_obj,
                        n_pcs = 20,
                        res = c(
                          0.005, 0.01, 0.04, 0.05,
                          seq(from = 0.1, to = 1, by = 0.1)
                        ),
                        select_genes = rownames(seur_obj),
                        elbow_dims = 50,
                        tsne = TRUE,
                        dir_lab = "all_celltypes",
                        save_dir = getwd(),
                        plotheight = 5,
                        plotwidth = 5,
                        sketch = FALSE,
                        sketch_ncells = 50000,
                        sketch_method = c("LeverageScore", "Uniform"),
                        sketched_assay_name =  "sketch") {
  out_dat_dir <- paste0(save_dir, "/data/processed_seurat/")
  if(!dir.exists(out_dat_dir)){
    dir.create(out_dat_dir, recursive = TRUE)
  }
  print("Normalising data.")
  seur_obj <- Seurat::NormalizeData(seur_obj)
  print("Looking for variable genes.")
  seur_obj <- Seurat::FindVariableFeatures(seur_obj)
  saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
  print("Scaling data.")
  scale_genes <- select_genes
  seur_obj <- Seurat::ScaleData(seur_obj, features = scale_genes)
  saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
  
  if(sketch == TRUE){
    print("Sketching data.")
    seur_obj <- SketchData(
      object = seur_obj,
      ncells = sketch_ncells,
      method = sketch_method,
      sketched.assay = sketched_assay_name
    )
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
    DefaultAssay(seur_obj) <- "sketch"
    # perform clustering workflow
    print("Find variable features of sketched data.")
    seur_obj <- Seurat::FindVariableFeatures(seur_obj)
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
    print("Scale sketched data.")
    seur_obj <- Seurat::ScaleData(seur_obj)
    print("Run PCA on sketched data.")
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
    seur_obj <- Seurat::RunPCA(seur_obj, assay = "sketch", reduction.name = "pca.sketch")
    print("Find nearest neighbors using sketched data.")
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
    seur_obj <- Seurat::FindNeighbors(seur_obj, assay = "sketch", reduction = "pca.sketch", 
                              dims = 1:n_pcs)
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
    print("Find clusters using sketched data.")
    seur_obj <- Seurat::FindClusters(seur_obj, cluster.name = "seurat_cluster.sketched", 
                             resolution = res)
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
    print("Run UMAP using sketched data.")
    seur_obj <- Seurat::RunUMAP(seur_obj, reduction = "pca.sketch", 
                       reduction.name = "umap.sketch", return.model = T, 
                       dims = 1:n_pcs)
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
  }else{
    print("Run PCA.")
    seur_obj <- Seurat::RunPCA(seur_obj, features = VariableFeatures(object = seur_obj))
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
  }
    
  e_plot <- Seurat::ElbowPlot(seur_obj, ndims = elbow_dims) +
      geom_vline(xintercept =  n_pcs,
                 color = "red", size=0.5)
  print(e_plot)
  save_plot_path <- paste0(save_dir, "/outs/", dir_lab, "/plots/elbow_plot")
  
  if(dir.exists(save_plot_path) == FALSE){
    dir.create(save_plot_path, recursive = TRUE)
  }
  
  pdf(paste0(save_plot_path, "/elbow.pdf"),
      height=plotheight, width = plotwidth)
  print(e_plot)
  dev.off()
  
  if(sketch == "FALSE"){
    print("Find nearest neighbours.")
    seur_obj <- Seurat::FindNeighbors(seur_obj, dims = 1:n_pcs)
    print("Find clusters.")
    seur_obj <- Seurat::FindClusters(seur_obj, resolution = res)
    print("Run UMAP.")
    seur_obj <- Seurat::RunUMAP(seur_obj, dims = 1:n_pcs)
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
    if (tsne == TRUE) {
      print("Run TSNE")
      seur_obj <- Seurat::RunTSNE(seur_obj, dims = 1:n_pcs)
    }
    saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
  }
  saveRDS(seur_obj, paste0(out_dat_dir, "processed_srt.RDS"))
  print("DONE")
  return(seur_obj)
  
}




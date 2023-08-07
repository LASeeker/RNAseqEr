#' Fast Seurat Processing of RNAseq data until dimensional reduction
#' @description This function performs the standard Seurat processing on a quality
#' controlled Seurat object an returns a Seruat object with added dimensional reduction
#' and clustering at multiple resolutions.
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
#' @return returns a Seruat object that is normalised, has information on linear and
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
                        tsne = TRUE) {
  seur_obj <- Seurat::NormalizeData(seur_obj)
  seur_obj <- Seurat::FindVariableFeatures(seur_obj)
  scale_genes <- select_genes
  seur_obj <- Seurat::ScaleData(seur_obj, features = scale_genes)
  seur_obj <- Seurat::RunPCA(seur_obj, features = VariableFeatures(object = seur_obj))
  print(Seurat::ElbowPlot(seur_obj, ndims = elbow_dims))
  seur_obj <- Seurat::FindNeighbors(seur_obj, dims = 1:n_pcs)
  seur_obj <- Seurat::FindClusters(seur_obj, resolution = res)
  seur_obj <- Seurat::RunUMAP(seur_obj, dims = 1:n_pcs)
  if (tsne == TRUE) {
    seur_obj <- Seurat::RunTSNE(seur_obj, dims = 1:n_pcs)
  }
  return(seur_obj)
}

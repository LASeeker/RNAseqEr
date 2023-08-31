#' Example dataset of oligodendrocyte nuclei derived from the human CNS
#'
#' @docType data
#'
#' @usage data(oligodendrocytes)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references Seeker, L.A., Bestard-Cuche, N., Jäkel, S. et al.
#' Brain matters: unveiling the distinct contributions of region, age, and
#' sex to glia diversity and CNS function. acta neuropathol commun 11, 84 (2023).
#' https://doi.org/10.1186/s40478-023-01568-z
#'
#'
#' @examples
#' library(Seurat)
#' data(oligodendrocytes)
#' DimPlot(oligodendrocytes, group.by = "cluster_id", cols = colour_palette())
"oligodendrocytes"

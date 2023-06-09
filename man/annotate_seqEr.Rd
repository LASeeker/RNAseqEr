% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_seqEr.R
\name{annotate_seqEr}
\alias{annotate_seqEr}
\title{Automatic annotation of broad celltypes using ScType}
\usage{
annotate_seqEr(
  seur_obj,
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
  dge = FALSE
)
}
\arguments{
\item{seur_obj}{Seurat object}

\item{ident}{meta data column that contains clustering resolution of interest
for broad celltype annotation}

\item{tissue_ref}{Specify which tissue is being used by chosing from the following:
"Brain", "Immune system", "Pancreas", "Liver", "Eye", "Kidney", "Brain", "Lung",
"Adrenal", "Heart", "Intestine", "Muscle", "Placenta", "Spleen", "Stomach",
"Thymus"}

\item{test_use}{statistical test used for differential gene expression analysis.
Default is "MAST".}

\item{min_pct}{minimum percentage of cells within the cluster of interest
expressing the gene. Default is 0.25.}

\item{logfc_threshold}{Minimal log2 fold change for differential gene expression.
Default is 0.8.}

\item{assay_use}{Which assay should be used for differential gene expression analusis.
Default is "RNA".}

\item{save_dir}{Directory to be used for saving output.}

\item{custom_ref_genes}{TRUE/FALSE whether a list of genes and celltypes is
being provided by the user.}

\item{cusom_gene_list}{if custom_ref_genes == TRUE provide custom gene list here.}

\item{customclassif}{which name should be used for the meta data column saving the
scType annotation information. Default is ScType_annotation"}

\item{plot_reduction}{which dimensional reduction should be used to plot data
with ScType labels. Default is "umap"}

\item{plot_label}{TRUE/FALSE whether ScType labels should be displayed on the
DimPlot. Default is TRUE.}

\item{plot_repel}{TRUE/FALSE whether ScType labels should be repelled on the
DimPlot. Default is TRUE.}

\item{plot_width}{width of DimPlot in output file. Default = 8.}

\item{plot_height}{height of DimPlot in output file. Default = 8.}

\item{dge}{TURE/FALSE whether differential gene expression should be perfomed.
Default is FALSE, because it may take some time.}
}
\value{
seurat object with added ScType annotation
}
\description{
ScType (https://www.nature.com/articles/s41467-022-28803-w) is used to identify
cell types automatically. If a onsiderable number of cells remain unknown,
an additional step which is based on differentially expressed genes is used
to attempt an celltype identification.
}
\examples{
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
pbmc_small <- annotate_seqEr(pbmc_small,
    tissue_ref= "Immune system",
    ident = "RNA_snn_res.1")
}

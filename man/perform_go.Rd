% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform_go.R
\name{perform_go}
\alias{perform_go}
\title{perfome_go}
\usage{
perform_go(
  seur_obj,
  gene_list,
  min_log2FC = 0.25,
  reverse = FALSE,
  translate_gene_id_from = "SYMBOL",
  translate_gene_id_to = "ENTREZID",
  org_use = "org.Hs.eg.db",
  ontology = "BP",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  read_able = TRUE,
  logFC_column = "avg_log2FC"
)
}
\arguments{
\item{seur_obj}{Seurat object}

\item{gene_list}{list of genes to use for gene ontology analyses}

\item{min_log2FC}{log2 fold cutoff for genes to be considered in gene ontology
analysis}

\item{reverse}{if genes smaller than the log2 FC cuttoff should be considered
instead, for example for smaler than -0.25}

\item{translate_gene_id_from}{current format of gene ids, default is "SYMBOL"}

\item{translate_gene_id_to}{specify which gene id's should be used (normally
ENTREZID), default is c("ENTREZID", "PATH", "GO", "ALIAS", "GENENAME")}

\item{org_use}{specify reference for your organism, default is human ("org.Dr.eg.db")}

\item{ontology}{what kind of gene ontology should be performed. Default is "BP" for
biological pathway but "MF" for molecular function for example is valid, too.}

\item{pvalue_cutoff}{p value cutoff for significant gene ontology results}

\item{qvalue_cutoff}{q value cutoff for significant gene ontology results}

\item{read_able}{gene names in GO output readable? Default is TRUE}
}
\value{
returns a gene ontology result (add data format) ready for plotting
and saving
}
\description{
Run gene ontology using cluster profiler
}
\examples{
data(dge_1_0)

perform_go(pbmc_small, gene_list = dge_1_0)

}

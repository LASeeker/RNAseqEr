% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_mark_list.R
\name{gen_mark_list}
\alias{gen_mark_list}
\title{Generate marker list}
\usage{
gen_mark_list(
  file_dir = getwd(),
  ad_pval = 0.05,
  avg_log = 1.2,
  pct_1 = 0.25,
  pct_2 = 0.6,
  pairwise = FALSE,
  n_top = 10,
  condition = FALSE,
  test_cond = NULL
)
}
\arguments{
\item{file_dir}{directory where differential gene expression results are saved
as .csv files.}

\item{ad_pval}{threshold for adjusted p-value, default is 0.05}

\item{avg_log}{threshold for minimum average log 2 fold change, default is 1.2}

\item{pct_1}{threshold for percentage of cells/ nuclei within the cluster
of interest expressing the gene. Default is 0.25}

\item{pct_2}{threshold for percentage of cells/ nuclei outside the cluster
of interest expressing the gene. Default is 0.6}

\item{pairwise}{TRUE/FALSE whether differential gene expression was performed
pairwise or not}

\item{n_top}{number of top genes per cluster considered. Default is 10.}

\item{condition}{TRUE/FALSE whether cluster marker or condition markers are
tested. Default is FALSE, so it is assumed that cluster markers are being
tested. TRUE enables the test of different conditions such as age, sex,
condition and treatment group separately.}

\item{test_cond}{FALSE by default possible to set to conditions of interest
that for which marker gene lists should be compiled. The conditions have to
appear in the Differential gene expression results that are used here as input
for example c("AgeGroup", "Tissue")}
}
\value{
filtered data frame of potential cluster marker genes with n_top candidate
genes per cluster/ group of interest
}
\description{
Use diffential gene expression analysis results to generate a list of top markers
for each cluster/ condition.
}
\examples{
library(Seurat)
library(dplyr)
library(here)
save_dir_cr <- here("data", "dge")
dir.create(save_dir_cr, recursive = TRUE)
df_0_8 <- as.data.frame(dge_0_8)
df_1_0 <- as.data.frame(dge_1_0)
write.csv(df_0_8, here(save_dir_cr, "dge_0_8.csv"))
write.csv(df_1_0, here(save_dir_cr, "dge_1_0.csv"))
clu_mark <- gen_mark_list(file_dir = here("data", "dge"))
}

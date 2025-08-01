# RNAseqEr

Fast and reproducible single cell/ nucleus RNAseq analysis

[![R-CMD-check](https://github.com/LASeeker/RNAseqEr/workflows/R-CMD-check/badge.svg)](https://github.com/LASeeker/RNAseqEr/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This package was written to allow for a fast and reproducible processing
of single cell or single nucleus RNA sequencing (sc/snRNAseq) data in R
by summarizing common workflows into fewer functions and by providing
solutions for bottle-neck steps that often slow down progress such as
the cell type annotation and the decision on a suitable clustering
resolution that is biologically meaningful.

The workflow starts with a pre-processed dataset where the following is
considered:
* genes that are not expressed may be removed from the dataset (gene QC)
* low quality cells are removed from the dataset
* ambient RNA may be removed
* doublets may be removed or flagged
* dataset may be batch corrected and/or integrated if required.

## Workflow

The workflow consists of:

1. **Performing standard Seurat processing** including:
   * Normalisation
   * finding of variable genes
   * linear and non-linear dimensional reductions
   * Clustering at different resolutions

2. **Finding a stable cluster resolution** for automatic celltype
   annotation using an adapted ScType algorithm

3. **Subsetting the datasets** for cell lineages

4. **Finding a suitable clustering resolution** for cell lineage datasets
   (bottleneck)

5. **Cluster quality control**

6. **Identification of best cluster marker genes** for validation

7. **Differential abundance** with variables of interest

8. **Differential gene expression** with variables of interest and
   identification of condition markers

9. **Gene ontology** based on cluster and condition markers

10. **Building of a shiny app**

11. **Writing of a summary report** that describes what has been done, what
    the results are and where to find them.

The processed Seurat objects are saved which enables the addition of
other tailored analyses.

All results should be thoroughly checked and manually corrected if
necessary.

## Installation

```r
# Install from GitHub
devtools::install_github("LASeeker/RNAseqEr")
```

## Quick Start

```r
# Load libraries
library(RNAseqEr)
library(here)
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(MAST)
library(dplyr)

# Load example data
data(cns)

# Multi-species support
result <- quick_RNAseqEr(mouse_data, species = "Mouse", tissue_ref_annotation = "Liver")
result <- quick_RNAseqEr(zebrafish_data, species = "Zebrafish", tissue_ref_annotation = "Brain")
```

## API Key Setup for LLM Integration

RNAseqEr supports LLM integration for manuscript generation. To use this feature, you need to set up your API keys.

### **Option 1: Using dotenv package (Recommended)**

1. **Install the dotenv package:**
   ```r
   install.packages("dotenv")
   ```

2. **Create a `.env` file** in your project root:
   ```
   API_KEY_GPT4=your-openai-api-key-here
   ANTHROPIC_API_KEY=your-anthropic-api-key-here
   ```

3. **Test your setup:**
   ```r
   library(RNAseqEr)
   setup_rnaseqer_env()
   ```

### **Option 2: Manual setup**

If you prefer to set the API key manually:

```r
library(RNAseqEr)
set_api_key_manual("your-api-key-here", "openai")
test_api_key_setup()
```

### **Usage in Analysis:**

```r
# Enable manuscript generation
result <- quick_RNAseqEr(data,
                         manuscript_generation = TRUE,
                         manuscript_llm_provider = "openai")
```

## Multi-Species Support

RNAseqEr supports analysis of single-cell RNA-seq data from multiple species:

### **Supported Species:**
- **Mammals:** Human, Mouse, Rat, Cow, Pig, Dog, Cat, Horse, Sheep, Goat, Rabbit, Guinea Pig, Hamster
- **Primates:** Monkey, Chimpanzee, Gorilla, Orangutan, Marmoset, Squirrel Monkey, Tarsier, Bushbaby
- **Other:** Zebrafish, Drosophila, Chicken

### **Automatic Adaptations:**
- **Gene Ontology Databases:** Automatically selects appropriate organism database
- **Tissue References:** Adapts ScType annotation for different tissues
- **Gene ID Mapping:** Handles species-specific gene identifiers
- **Statistical Methods:** All statistical methods work across species

### **Usage Examples:**
```r
# Human brain analysis (default)
result <- quick_RNAseqEr(human_data, species = "Human", tissue_ref_annotation = "Brain")

# Mouse liver analysis
result <- quick_RNAseqEr(mouse_data, species = "Mouse", tissue_ref_annotation = "Liver")

# Zebrafish brain analysis
result <- quick_RNAseqEr(zebrafish_data, species = "Zebrafish", tissue_ref_annotation = "Brain")

# Drosophila analysis
result <- quick_RNAseqEr(drosophila_data, species = "Drosophila", tissue_ref_annotation = "Brain")
```

## Test with provided test dataset or own data

The package comes with a dataset from the human CNS (Seeker et al. 2023)
that includes different cell types from different CNS regions and
different donor age and sex groups. This data can be loaded with:

```r
cns
```

## Standard Seurat processing

It can be seen that the cns data above already contains normalized and
scaled data as well as different dimensional reductions. Your data will
not have any of those at this point of the analysis.
[Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) is
a fantastic and very user-friendly tool for adding those information to
your data and I highly recommend reading through their papers and
vignettes. However, if you analyse many datasets, you will also see that
your work is going to be very repetitive, using the same Seurat
functions again and again without tweaking them necessarily a lot. So
here we present a function that automatically runs Seurat's standard
normalisation, selection of variable genes, dimensional reduction and
clustering at different resolutions.

```r
cns <- seurat_proc(cns)
```

The function plots an elbow plot so that you can check if the standard of 20 principle
components is suitable for your data. If you think you need to use less or more 
(20 is far away from the elbow in the plot), you can re-run the function while
providing an additional argument:

```r
cns <- seurat_proc(cns, n_pcs = 25)
```

## Plotting of different clustering resolutions

In the previous step you created different clustering resolutions, but
which is the most appropriate one for your data? This is often a tricky
question to answer and the literature on this is vast.

A first good step is to visualise different clustering resolutions. You
can easily do this with the function below which will also save
dimensional reduction plots (default is umap, tsne can be used with
use_reduction = "umap") to a folder. Remember to set your save_directory
that you'd like to use, the default is the current working directory.

```r
plot_list(cns, save_dir = "../")
```

## Calculate cluster purity at different resolutions

### RNAseqEr approach based on purity measures

The intention above was to create clustering from a too low to a too high resolution. 
But which one is appropriate? This depends in this work for on if you have a 
very homogeneous cell population or a mix of different cell types (which is the
case here), where we would like to do first a rough annotation.

If the second is true, we will use cluster purity measures to first separate the
dataset best for clusters that are physically distinct on in the umap/tsne space 
and thus most likely represent different cell types.

```r
clu_pure(cns, save_dir = "../")
```

```r
keep_res <- max_pure(save_dir = "../", dir_lab = "all_celltypes")

DimPlot(cns, group.by = keep_res, cols = colour_palette())
```

### Using CHOIR
A recent publication by Sant et al. proposes 

Sant, C., Mucke, L. & Corces, M.R. CHOIR improves significance-based detection of cell types and states from single-cell data. Nat Genet (2025). https://doi.org/10.1038/s41588-025-02148-8

```r
# Note: CHOIR is not included in the package dependencies
# cns_choir <- CHOIR(cns)
# cns_choir <- runCHOIRumap(cns_choir, reduction = "P0_reduction")
# plotCHOIR(cns_choir)    
# plotCHOIR(cns_choir, accuracy_scores = TRUE, plot_nearest = FALSE)
```

## Cell type annotation

Above we found a cluster resolution that separates what could be
celltypes, however, we don't know which cell types are present. We are
therefore using an annotation based on
[ScType](https://www.nature.com/articles/s41467-022-28803-w). We are
looking for a rough annotation, that can be used for subsetting the
dataset into lineages for subsequent further clustering. SCType
sometimes returns "Unknown" and we therefore adapted it for those cases
to not only consider gene expression but also differential gene
expression. Running a differential gene expression analysis at this step
is time consuming but also valuable as it will help to double check the
automatic annotation. The differential gene expression results are saved
to file and can be used for plotting heatmaps for example that show cell
type segregation.

N.B. The SCType and RNAseqER annotation will look identical as long as
there are no "Unknown" in the former.

```r
cns <- annotate_seqEr(cns, ident_use = keep_res, tissue_ref = "Brain",
                      save_dir = "../")

cns <- annotate_seqEr(cns, ident_use = "RNA_snn_res.0.8", tissue_ref = "Brain",
                      save_dir = "../")
```

The annotation above looks very close to the manual annotation we used in 
Seeker et al. 2023:

```r
DimPlot(cns, group.by = "rough_annot", cols = colour_palette(), label = TRUE)
```

## Finer annotation

We may or may not be interested in all celltypes that a dataset
contains. In many cases it is beneficial to subset the dataset into cell
lineage datasets that then can be used for finer clustering and
downstream analyses. The following function subsets the dataset for all
clusters that have been identified and annotated above using the
RNAseqEr annotation. Any metadata column can be used for subsetting.
Each cell lineage will be saved in a separate .RDS file.

```r
subset_RNAseqEr(cns, save_dir = here())

list.files(here("outs", "data"))
```

## Finer clustering

Finding an appropriate clustering resolution within a cell lineage is
often a step that takes a long time and is associated with insecurity,
particularly when the data scientist is yet inexperienced. From a
biological point of view we seek to find a resolution with the maximum
number of clusters that still allows the identification of cluster
markers that can be used for validation using other lab techniques such
as immuno-fluorescence or RNAscope. We show here our approach using the
example of oligodendrocytes. We also provide a function using the same
approach for all the subsetted datasets automatically.

```r
seur <- readRDS(here("outs", "data", "Oligodendrocytes.RDS"))
seur <- seurat_proc(seur, tsne = FALSE)
```

```r
plot_list(seur, save_dir = "../", dir_lab = "Oligodendrocytes")
```

In the above plots it can be seen that some resolutions do not seem to
make sense as they either don't cluster at all (< 0.1), have the same
clusters as another resolution (for example 0.6 and 0.7) or simply seem
to overcluster, where cluster separations are visually not detectable
and cluster labels overlap completely (>1).

All other cluster resolutions may make sense, so how to pick a suitable
one? Our approach is to decide by looking at the results of
differential gene expression analyses to see whether or not we can
identify markers for most clusters. Why most and not all? We noticed
that in many datasets, one cluster seems to be the "average" that does
express all the cell lineage markers but nothing special over and above.
It may be reasonable to keep a clustering resolution where a cluster is
negative for markers that characterise other clusters.

So all potentially interesting resolutions can be used below:

```r
# Pick interesting resolutions
int_resol <- c(0.1, 0.2, 0.3, 1, 1.3)
int_res <- paste0("RNA_snn_res.", int_resol)
```

```r
all_res_mark <- int_res_all_mark(seur_obj = seur,
                                 int_cols = int_res,
                                 save_dir = save_dir,
                                 dir_lab = "Oligodendrocytes")
```

```r
# see the output 
list.files(here("outs", "Oligodendrocytes", "tables", "cluster_marker", "overall"))
```

The function above generates two outputs for each resolution, one set
for unfiltered ("all") and one for further filtered ("fil") values that
can be tweaked in the function arguments.

We also add a pairwise comparison of clusters:

```r
pairwise_dge(seur_obj = seur,
             int_cols = int_res,
             save_dir = save_dir,
             dir_lab = "Oligodendrocytes")
list.files(here("outs", "Oligodendrocytes", "tables", "cluster_marker", "pairwise"))
```

Above we generate a whole lot of output that nobody wants to comb through 
manually. So below we are going to read in all the results and summarise them
in a meaningful way. 

The two functions below look into the results of the overall (first) and pairwise
(second) differential gene expression analysis results and keeps the top n 
markers (default n = 10) per cluster. 

Those genes are then concatenated and duplicates are filtered out, leaving the
user with a condensed list of the most interesting potential cluster markers which
are used to determine the similarity/difference of clusters to find an appropriate
clustering resolution.

```r
clu_mark <- gen_mark_list(file_dir = here("outs", 
                                          "Oligodendrocytes", 
                                          "tables", 
                                          "cluster_marker", 
                                          "overall"))

clu_mark_pw <- gen_mark_list(file_dir = here("outs", 
                                          "Oligodendrocytes", 
                                          "tables", 
                                          "cluster_marker", 
                                          "pairwise"),
                             pairwise = TRUE)

int_genes <- c(clu_mark$gene, clu_mark_pw$gene)
int_genes <- unique(int_genes)
```

## Identifying clusters that are too similar to find most appropriate resolution

Below we use heatmaps to visually screen differences and similarities
between clusters. The function also uses average expression values to
compare each cluster pair at different resolutions pairwise. If they are
too similar based on their maximum difference, mean difference or
Euclidean distance, it's recommended to merge clusters. Each of those
incidences adds 1 to a cluster similarity score. If a resolution has a
high overall cluster similarity score it means that resolution includes
a lot of clusters that are very similar. The cluster resolution with the
lowest cluster similarity score is the one the user should consider.

```r
hm_dir <- here("outs", "Oligodendrocytes", "plots", "heatmaps")

dir.create(hm_dir, recursive = TRUE)

df_seqEr <- heatmap_seqEr(seur,
                          use_resol = FALSE,
                          col_names = int_res,
                          int_genes = int_genes,
                          save_dir = hm_dir,
                          max_diff_threshold = 10,
                          mean_diff_thres = 0.1)
```

To look at the cluster similarity scores, inspect the dataframe below. It shows that
in this case the resolution 1 may be the most appropriate. More generally we use the 
highest resolution with the lowest cluster similarity score (tail(x,1)).

```r
df_seqEr
```

The below only works if there was some over-clustering and there is
indeed one cluster with a minimum cluster similarity score.
Alternatively this can be determined manually.

```r
# res_df <- subset(df_seqEr, df_seqEr$cluster_sim_score == min(df_seqEr$cluster_sim_score))
# keep_res <- tail(res_df$resolution,1)
```

```r
keep_res <- "RNA_snn_res.1"
```

```r
Idents(seur) <- keep_res

seur$cluster_id <- paste0("Oligodendrocytes_", seur@meta.data[[keep_res]])

DimPlot(seur, group.by = "cluster_id", label = TRUE, cols = colour_palette())
```

## Cluster Quality control

Sometimes individual samples can cluster separately from the rest which is 
more likely to be due to technical reasons than showing a biologically interesting
signal. This may particularly happen when working with human or other non-model
organism data where the post-mortem conditions are more variable and there is more
genetic variation and variation in lifestyle, disease burden and medication. 

We always perform some form of cluster quality control to capture clusters that
are made up only by a few samples or donors. If tiny clusters are as a consequence 
removed from the dataset, we keep the rest of the dataset as it is. If a significant 
proportion is removed, it is a safer option to then repeat the selection of variable 
genes up to the re-clustering. 

The function below can be used on any metadata column to test for cluster
biased clustering. 

```r
seur <- cluster_qc(seur_obj = seur,
                   save_dir = save_dir,
                   cluster_col = "cluster_id",
                   vars = c("process_number", "caseNO", "Tissue"),
                   dir_lab = "Oligodendrocytes")
```

## Differential abundance testing

Often it is of interest to see if the number of cell/nuclei per cluster
vary depending on the condition for example. Different methods have been
implemented to test this question in a way that is statistically sound
which is less trivial than it may appear at first. In our hands Milo
works well and we implemented it here. We encourage the user to try
and test other methods as well.

```r
milo_obj <- abundance_test(seur_obj = seur,
                           cluster_col = "cluster_id",
                           sample_id = "uniq_id",
                           cols_interest = c("uniq_id",
                                             #"Tissue",
                                             #"gender",
                                             "AgeGroup" #,
                                             #"caseNO"
                                             ),
                           test_factors = c(#"Tissue", 
                                            #"gender", 
                                            "AgeGroup"),
                           dir_lab = "Oligodendrocytes",
                           k = 10,           # Reduced from 30
                           d = 10,           # Reduced from 30
                           prop = 0.1,       # Reduced from 0.2
                           refined = FALSE  #
                           )
```

Differential abundance testing makes more sense for larger datasets. Therefore,
I am showing it below for the cns dataset containing all cell types.

```r
milo_obj <- abundance_test(seur_obj = cns,
                       cluster_col = "rough_annot",
                       sample_id = "uniq_id",
                       cols_interest = c("uniq_id",
                                          #"Tissue",
                                          #"gender",
                                          "AgeGroup"#,
                                          #"caseNO"
                                          ),
                      test_factors = c(#"Tissue", 
                                       #"gender", 
                                       "AgeGroup"
                                       ),
                      save_dir = "../",
                           k = 10,           # Reduced from 30
                           d = 10,           # Reduced from 30
                           prop = 0.1,       # Reduced from 0.2
                           refined = FALSE  #
                      )
```

## Differential gene expression with condition

Next it is interesting to find out which genes are differentially
expressed with variables of interest such as disease condition,
treatment group, age group or tissue region for example. One approach is
to see which genes are expressed in the clusters more abundant in those
groups (see above). The other is to look within cell lineage or
clusterwise for genes that are differentially expressed.

```r
RNAseqEr::find_cond_markers(seur,
                  int_cols = c("AgeGroup","Tissue", "gender"),
                  cluster_id = "cluster_id",
                  dir_lab = "Oligodendrocytes",
                  save_dir = "../")
            
```

The subsetted oligodendrocyte dataset is very small to test for
differentially expressed genes. This is why we show the workflow
downstream with the cns dataset.

```r
RNAseqEr::find_cond_markers(cns,
                  int_cols = c("AgeGroup","Tissue", "gender"),
                  cluster_id = "rough_annot",
                  save_dir = "../")
            
```

The example below shows how to collect data from the cluster-wise
differential gene expression analysis that we performed above and plot
them as violin plots.

```r
library(here)
cond_genes <- gen_mark_list(file_dir = here("outs", 
                                            "all_celltypes", 
                                            "tables", 
                                            "condition_mark" , 
                                            "rough_annot" , 
                                            "clusterwise"),
                            condition = TRUE,
                            test_cond = c("AgeGroup", "Tissue", "gender"))

# Subset for the condition we would like to plot starting with age
age_mark <- subset(cond_genes, 
                              cond_genes$condition == "AgeGroup")

#There is no point in plotting the same genes several times. Therefore, we are
# removing duplicated gene names.
age_mark <- age_mark[!duplicated(age_mark$gene),]
head(age_mark)

# Save violin plots to file
save_vln(cns, 
         marker_list = age_mark$gene,
         save_label = "Age_mark",
         split_by = "AgeGroup",
         group_by = "rough_annot",
         plotheight = 50,
         save_dir = "../")
```

The same is also possible to perform on the overall differential gene expression
analysis from above by pointing at the corresponding folder.

```r
cond_genes_ov <- gen_mark_list(file_dir = here("outs", 
                                            "all_celltypes", 
                                            "tables", 
                                            "condition_mark" , 
                                            "rough_annot" , 
                                            "overall"),
                            condition = TRUE,
                            test_cond = c("AgeGroup", "Tissue"))

# Subset for the condition we would like to plot, this time choosing "Tissue"
tissue_mark <- subset(cond_genes_ov, 
                              cond_genes_ov$condition == "Tissue")

#There is no point in plotting the same genes several times. Therefore, we are
# removing duplicated gene names.
tissue_mark <- tissue_mark[!duplicated(tissue_mark$gene),]
head(tissue_mark)

# Save violin plots to file
save_vln(cns, 
         marker_list = tissue_mark$gene,
         save_label = "Tissue_mark",
         split_by = "Tissue",
         group_by = "rough_annot",
         plotheight = 30,
         save_dir = "../")
```

Similarly, Feature plots can be useful to visualize differences in
gene expression. We can use the gene list above to visualize differences with a
condition group such as tissue.

```r
save_feat_plots(cns, 
                marker_list = tissue_mark$gene,
                save_label = "Tissue_mark",
                split_by = "Tissue",
                plotheight = 30,
                save_dir = "../")
```

## Gene ontology analysis

It may also be useful to check whether there are certain pathways differentially 
regulated in certain clusters or with a certain condition. There are many 
tools available for gene ontology or gene set enrichment analyses and even more
ways to visualize the results. We use ClusterProfiler below and generate an output 
that can be easily plotted and saved as wished. 

Please note that the usage below is just for the purpose of demonstrating the 
function. When you use the function for your purposes, make sure to filter 
your list of genes so that you retain only those that are significantly upregulated
in your group of interest. Run the function for each condition/cluster. 

```r
go_results <- perform_go(cns,
                         gene_list = tissue_mark,
                         min_log2FC = 0.25,
                         reverse = FALSE,
                         translate_gene_id_from = "SYMBOL",
                         translate_gene_id_to = "ENTREZID",
                         org_use = "org.Hs.eg.db",
                         ontology = "BP",
                         pvalue_cutoff = 0.05,
                         qvalue_cutoff = 0.05,
                         read_able = TRUE)

dotplot(go_results)
```

## Output Structure

RNAseqEr creates a comprehensive output structure organized in the `outs/` directory. Here's the complete folder structure:

```
outs/
├── all_celltypes/                    # Main analysis results
│   ├── plots/
│   │   ├── resolution_plots/         # UMAP plots at different resolutions
│   │   ├── final_clustering/         # Final clustering plots
│   │   ├── Condition/                # Condition-specific plots
│   │   │   ├── AgeGroup_mark/
│   │   │   │   └── violin_plots/
│   │   │   └── Tissue_mark/
│   │   │       └── violin_plots/
│   │   └── heatmaps/                 # Cluster similarity heatmaps
│   ├── tables/
│   │   ├── DGE/
│   │   │   └── broad_celltype_markers/  # Differential gene expression results
│   │   ├── condition_mark/
│   │   │   └── RNAseqEr_annotation/
│   │   │       ├── clusterwise/      # Cluster-wise condition markers
│   │   │       └── overall/          # Overall condition markers
│   │   └── cluster_marker/
│   │       ├── overall/              # Overall cluster markers
│   │       └── pairwise/             # Pairwise cluster comparisons
│   └── summary_figures/
│       └── conditions/
│           └── volcano/               # Volcano plots for condition markers
├── data/                             # Processed Seurat objects
│   ├── processed/
│   │   └── all_data/                 # Main processed dataset
│   ├── Astrocytes.RDS                # Subsetted cell lineage data
│   ├── Oligodendrocytes.RDS
│   ├── Neurons.RDS
│   └── ...
├── supplementary_tables/              # Summary tables
│   ├── overall_condition_markers.csv
│   └── clusterwise_condition_markers.csv
├── shiny_app/                        # Interactive Shiny application
│   ├── ui.R
│   ├── server.R
│   └── www/
├── reports/                          # Analysis reports
│   ├── methods_report.md
│   ├── decision_justifications.md
│   └── complete_analysis_report.md
└── manuscript/                       # LLM-generated manuscript (optional)
    ├── manuscript_draft.md
    ├── figure_selection.md
    └── summary_statistics.md
```

### Key Output Files:

**Plots:**
- `resolution_plots/`: UMAP plots at different clustering resolutions
- `final_clustering/`: Final clustering plots for each cell lineage
- `Condition/`: Violin plots and feature plots for condition markers
- `heatmaps/`: Cluster similarity heatmaps
- `volcano/`: Volcano plots showing differential expression

**Tables:**
- `DGE/`: Differential gene expression results
- `condition_mark/`: Condition-specific marker genes
- `cluster_marker/`: Cluster-specific marker genes
- `supplementary_tables/`: Summary CSV files

**Data:**
- `.RDS` files: Processed Seurat objects for each cell lineage
- `processed/`: Main processed dataset

**Applications:**
- `shiny_app/`: Interactive web application for data exploration
- `reports/`: Comprehensive analysis reports
- `manuscript/`: LLM-generated manuscript draft and figure selection

## Generation of a Shiny app

```r
create_shiny(cns,
             shiny_name = "all_celltypes",
             save_dir = "../",
             read_file = FALSE,
             default_1 = "AgeGroup",
             default_2 = "Tissue",
             assay_use = "RNA",
             gex_slot = c("data", "scale.data", "counts"),
             gene_mapping = FALSE,
             default_gene1 = "MALAT1",
             default_gene2 = "GAPDH",
             default_multigene = NA,
             default_dimred = c("umap1", "umap2"),
             author = "Luise A. Seeker",
             title = "RNAseqEr: Fast and reproducible sc/snRNAseq data analysis",
             journal = "Bioinformatics",
             volume  = "1",
             page    = "1",
             year    = "2024",
             doi     = "https://doi.org/10.1093/bioinformatics/xxx",
             link    = "https://github.com/LASeeker/RNAseqEr",
             shiny_title = "RNAseqEr Single Cell Analysis")
```

This function will generate a directory called "shiny_app" in the
provided save_dir. In that directory is a script called ui.R. When this
script is opened in R Studio, an option to "Run App" will appear on the
top right. When this is clicked, a website will be opened that can
published to share data with collaborators (secret link) and/or the
world. Several datasets can be visualized in the same shiny app as shown
in https://seeker-science.shinyapps.io/shiny_app_multi/ for example.
More data visualization modes are available in the drop down menu of
each dataset/cell type tab.

## Citation

If you use RNAseqEr in your research, please cite:

```r
citation("RNAseqEr")
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact

- **Author**: Luise A. Seeker
- **Email**: seeker.luise@gmail.com
- **GitHub**: [@LASeeker](https://github.com/LASeeker)
- **Issues**: [GitHub Issues](https://github.com/LASeeker/RNAseqEr/issues)

## Session Info

```r
sessionInfo()
```

# Comprehensive Seurat Single-Cell RNA-seq Pipeline

## Description

This script offers a robust and customizable command-line solution for processing single-nucleus RNA sequencing (snRNA-seq) data using the powerful Seurat R package. Designed for both raw 10X Genomics H5 data and existing Seurat objects, it guides users through essential analysis steps including rigorous quality control (QC), data normalization and scaling, dimensionality reduction (PCA), graph-based clustering, UMAP visualization, and optional differential gene expression (DEG) analysis for cluster marker identification. A key feature is the integrated support for doublet removal using output from the SouporCell tool. The pipeline is built for ease of use in computational environments and supports batch processing workflows.

## Prerequisites

This script requires R (4.3.1) and several R packages. Ensure R is installed on your system. The script includes a section at the beginning to check for and automatically install the necessary R packages if they are not already present. Ensure your R environment has access to the internet during the initial run if installations are needed.

The required packages are:
* `Seurat` (>= 5.2.1)
* `tidyverse` (>= 2.0.0)
* `patchwork` (>= 1.3.0)
* `argparse` (>= 2.2.5)
* `scales` (>= 1.3.0)
* `ggrepel` (>= 0.9.6)

## Installation

1.  Save the R script code (provided separately in a file named `seurat_pipeline.R`) into a text file.
2.  Make the script executable in your terminal:

    ```bash
    chmod +x seurat_pipeline.R
    ```

## Usage

The script is executed from the command line. Navigate to the directory where you saved the script in your terminal. Use the `-h` or `--help` flag to display a detailed help message with all arguments and their descriptions:

```bash
./seurat_pipeline.R -h
```

- --rds_input CHARACTER Path to an existing Seurat object saved as an RDS file. Use this if you are resuming a pipeline or starting with pre-processed data. (Either this or --h5_input is required)
- --h5_input CHARACTER Path to a raw 10X Genomics filtered feature-barcode matrix in H5 format (e.g., filtered_feature_bc_matrix.h5 from Cell Ranger). Use this for initial processing from raw data. (Either this or --h5_input is required)
- --project_name CHARACTER A name for your project. This name will be used in plot titles and as a prefix for all output files. (Default: SeuratProject)
- --soupor_cell_doublet_input CHARACTER Path to a SouporCell clusters.tsv file for doublet filtering. If provided, the script will filter out cells classified as 'doublet'. (Optional)
- --min_features INTEGER Minimum number of genes detected in a cell. Cells with fewer features will be removed during QC if --h5_input is used. (Default: 200)
- --min_counts INTEGER Minimum number of total UMI counts per cell. Cells with fewer counts will be removed during QC if --h5_input is used. (Default: 1000)
- --n_variable_features INTEGER The number of top variable features to identify and use for downstream steps like PCA and scaling. (Default: 2000)
- --normalization_method {LogNormalize,RC,CLR} Specifies the normalization method to apply. LogNormalize: Log-normalization by total library size. RC: Relative counts (counts divided by total counts). CLR: Centered log-ratio transformation. - (Default: CLR)
- --scaling_method {negbinom,linear} Specifies the statistical model used when scaling data. 'negbinom' is recommended for UMI data to better account for variance. (Default: negbinom)
- --pca_dimensions INTEGER The number of principal components to compute and retain for downstream analysis (e.g., FindNeighbors, RunUMAP). (Default: 15)
- --clustering_resolution DOUBLE A parameter influencing the granularity of clustering. Higher values typically result in more clusters. Experimentation is often needed to find an appropriate value for your dataset. (Default: 0.5)
- --clustering_algorithm {louvain,multilevel,leiden,slm} The algorithm used for graph-based clustering. louvain / multilevel: Standard Louvain algorithm. leiden: Leiden algorithm, generally recommended for larger datasets as it tends to produce more connected partitions. slm: SLM (Sequential Louvain Method). (Default: leiden)
- --working_dir CHARACTER The absolute or relative path to the directory where all output files will be saved. It is highly recommended to specify a dedicated directory for each run or sample. (Default: . - the current directory)
- --find_markers LOGICAL A boolean flag (TRUE or FALSE) to determine whether to run the FindAllMarkers step to identify genes differentially expressed in each cluster. Finding markers can be computationally intensive. (Default: TRUE)


## Example Commands
1. Running with a raw H5 file using default QC thresholds and all other default parameters:

```bash

./seurat_pipeline.R --h5_input /path/to/your/filtered_feature_bc_matrix.h5 --project_name MyProject
```

2. Running with a raw H5 file, specifying custom QC thresholds and a dedicated working directory:

```bash

./seurat_pipeline.R --h5_input /path/to/your/filtered_feature_bc_matrix.h5 --project_name MyProject --min_features 300 --min_counts 1500 --working_dir /path/to/output/directory
```

3. Running with a Seurat RDS file and disabling marker finding:

```bash

./seurat_pipeline.R --rds_input /path/to/your/seurat_object.rds --project_name ProcessedData --find_markers FALSE --working_dir /path/to/output/directory
```

4. Run multiple samples in one run:
```bash
#!/bin/bash

# Define the commands for each sample
commands=(
  "./seurat_pipeline.R --h5_input /path/to/sample1/h5 --soupor_cell_doublet_input /path/to/sample1/doublets.tsv --working_dir /path/to/output/sample1 --project_name Sample1 --find_markers FALSE"
  "./seurat_pipeline.R --h5_input /path/to/sample2/h5 --soupor_cell_doublet_input /path/to/sample2/doublets.tsv --working_dir /path/to/output/sample2 --project_name Sample2 --find_markers FALSE"
  "./seurat_pipeline.R --h5_input /path/to/sample3/h5 --soupor_cell_doublet_input /path/to/sample3/doublets.tsv --working_dir /path/to/output/sample3 --project_name Sample3 --find_markers FALSE"
  # Add commands for other samples
)

# Determine the number of parallel jobs (e.g., number of CPU cores)
num_jobs=$(nproc) # Or specify a number, e.g., num_jobs=8
echo "Running ${#commands[@]} commands using <span class="math-inline">num\_jobs parallel jobs\.\.\."
\# Execute the commands in parallel
printf "%s\\n" "</span>{commands[@]}" | parallel -j "$num_jobs"

echo "All parallel pipeline runs have finished."
```

## Pipeline Overview

### 1. Data Loading
- Load raw 10X data (`Read10X_h5`) or existing Seurat RDS (`readRDS`).
- Create Seurat object if starting from raw data.

### 2. Quality Control (QC)
- Compute mitochondrial percentage (`PercentageFeatureSet`).
- Filter cells by features/count thresholds.
- Optional doublet removal via SouporCell.
- Generate QC plots (Violin, FeatureScatter, Histograms).

### 3. Normalization & Scaling
- Normalize data (`NormalizeData`).
- Identify variable features (`FindVariableFeatures`).
- Scale data (`ScaleData`) using linear or negative binomial methods.

### 4. Dimensionality Reduction & Clustering
- Run PCA (`RunPCA`), plot results, and generate elbow plots.
- Construct cell-cell nearest neighbor graphs (`FindNeighbors`).
- Cluster cells (`FindClusters`) using chosen algorithm (e.g., Leiden).
- Visualize clusters with UMAP (`RunUMAP`).

### 5. Optional Marker Detection
- Find differentially expressed genes per cluster (`FindAllMarkers`).
- Default test: Wilcoxon; optionally use DESeq2 (suited for larger clusters).

### Final Output
- QC, PCA, UMAP plots, marker gene CSV.
- Processed Seurat object saved as `.rds`.

For detailed function documentation, see [Seurat v5 Docs](https://satijalab.org/seurat/).


## Output Files
All output files are saved in the directory specified by the --working_dir argument. The filenames are prefixed with the --project_name.

- [project_name]_QC_VlnPlots_PreFilter_Combined.png: Static combined VlnPlots showing distributions of nFeature_RNA, nCount_RNA, and percent.mt before any filtering. Includes summary statistics in subtitles.
[put name of QC_VlnPlots_PreFilter_Combined.png image here]

- [project_name]_QC_FeatureScatter_PreFilter_Combined.png: Static combined FeatureScatter plots showing nCount_RNA vs nFeature_RNA and nCount_RNA vs percent.mt before filtering, with log10 scale on the x-axis.
[put name of QC_FeatureScatter_PreFilter_Combined.png image here]

- [project_name]_QC_UMI_Histogram_PreFilter.png: Static histogram showing the distribution of UMI counts per cell before filtering. Includes summary statistics in the subtitle.
[put name of QC_UMI_Histogram_PreFilter.png image here]

- [project_name]_QC_VlnPlots_PostFilter_Combined.png: (Generated only if H5 input was used and cells remain after filtering) Static combined VlnPlots showing distributions after filtering. Includes summary statistics and filter information in subtitles.
[put name of QC_VlnPlots_PostFilter_Combined.png image here]

- [project_name]_QC_FeatureScatter_PostFilter_Combined.png: (Generated only if H5 input was used and cells remain after filtering) Static combined FeatureScatter plots after filtering, with log10 scale on the x-axis and filter information in subtitles.
[put name of QC_FeatureScatter_PostFilter_Combined.png image here]

- [project_name]_QC_UMI_Histogram_PostFilter.png: (Generated only if H5 input was used and cells remain after filtering) Static histogram showing the distribution of UMI counts per cell after filtering. Includes summary statistics and filter information in the subtitle.
[put name of QC_UMI_Histogram_PostFilter.png image here]

- [project_name]_QC_cells_summary.png: Bar plot visualizing the total number of cells before and after all filtering steps (SouporCell and/or min_features/counts if applicable).
[put name of QC_cells_summary.png image here]

- [project_name]_QC_reads_summary.png: Bar plot visualizing the total number of reads (UMI counts) before and after all filtering steps.
[put name of QC_reads_summary.png image here]

- [project_name]_Variable_Features.html: An interactive HTML plot of the identified highly variable features. You can open this file in a web browser to hover over points and see the gene names.
[put name of Variable_Features.html preview image here]

- [project_name]_PCA.1.png: Combined static plot showing the PCA cell embedding colored by cluster and an Elbow plot to assess PC significance.
[put name of PCA.1.png image here]

- [project_name]_PCA.2.png: Static heatmap showing the expression of the top genes contributing to the first few principal components.
[put name of PCA.2.png image here]

- [project_name]_neighbors_graph.rds: An R Data Serialization (RDS) file containing the igraph object representing the cell-cell nearest neighbor graph constructed by FindNeighbors().

- [project_name]_shared_nearest_neighbors_graph.rds: An R Data Serialization (RDS) file containing the igraph object representing the Shared Nearest Neighbor (SNN) graph.

- [project_name]_UMAP.png: Static UMAP plot visualizing the cell embedding in 2D, colored and labeled by the identified clusters.
[put name of UMAP.png image here]

- [project_name]_markers_res[resolution]_alg[algorithm].csv: (Generated only if --find_markers TRUE) A CSV file containing the table of differentially expressed genes for each cluster. The filename includes the clustering resolution and algorithm used.
[put name of markers.csv preview image here]

- [project_name]_processed.rds: The final processed Seurat object, saved as an RDS file. This object contains all the data and analysis results and can be loaded into R or RStudio for further interactive exploration.

### Doublet Filtering (SouporCell)
Use `--soupor_cell_doublet_input` to filter doublets identified by [SouporCell](https://github.com/wheaton5/souporcell). Retains cells labeled as 'singlet' or 'unclassified' and removes explicit 'doublets'.

### Optional Marker Finding
Activate differential expression analysis with `--find_markers TRUE`. Default is DESeq2 test.


### Working Directory
It is strongly recommended to use the --working_dir argument to specify a dedicated output directory for each sample or analysis run. This helps keep your project files organized and prevents output files from different runs from overwriting each other. The script will create the directory if it doesn't exist. Ensure the user running the script has write permissions to this directory.

## Tips and Best Practices
##### Choosing QC Thresholds: 
The default --min_features (200) and --min_counts (1000) are common starting points, but these should be adjusted based on the specific distribution of your data visible in the pre-filter QC plots and your knowledge of the expected cell quality. Look for inflection points in the cumulative distribution of counts and features.
#### Selecting PCA Dimensions: 
The --pca_dimensions argument is crucial for capturing the major sources of variation while reducing noise. The Elbow plot ([project_name]_PCA.1.png) can help you decide how many PCs to include. Look for an "elbow" point where the variance explained by subsequent PCs drops significantly. You can also consider using JackStraw analysis interactively in R/RStudio after the script runs to determine statistically significant PCs.
#### Choosing Clustering Resolution: 
The --clustering_resolution parameter directly impacts the number and granularity of the resulting clusters. Higher values lead to more clusters. Experiment with different resolutions to find a level of clustering that aligns with your biological questions and expected cell populations.
#### Interpreting Plots: 
Use the generated plots to assess the quality of your data and the results of each analysis step. The QC plots help you understand the initial data quality and the impact of filtering. The PCA, UMAP, and Variable Feature plots help you assess the dimensionality reduction and feature selection. The UMAP plot with clusters is the primary visualization for cell type identification.
#### Resource Management (Parallel Runs): 
When running the parallel Bash script for multiple samples, carefully consider the number of --pca_dimensions and whether --find_markers is enabled, as these steps can be computationally and memory intensive. Adjust the num_jobs parameter in your Bash script based on the available CPU cores and memory on your system to avoid overloading it.
#### Interactive Exploration: 
The final [project_name]_processed.rds Seurat object is your primary resource for further interactive analysis in R or RStudio. You can load this object to perform manual cluster annotation, visualize specific gene expression, run additional analyses, etc.

### Troubleshooting
##### Missing Packages: 
If the script fails with a missing package error, ensure your internet connection is active and R can access CRAN to install them. If behind a firewall, you might need to configure R's proxy settings.
##### File Not Found: 
Double-check the paths provided for --rds_input, --h5_input, and --soupor_cell_doublet_input. Ensure the files exist at those exact locations and that the user running the script has read permissions.
##### Seurat Function Errors: 
Errors during Seurat function calls (like NormalizeData, RunPCA, FindClusters) might indicate issues with the data itself (e.g., no cells remaining after filtering, all values being zero). Check the logging messages for specific error details and review the QC plots and summary statistics to understand the state of the data before the failing step.
##### Memory Errors: 
If the script runs out of memory, especially during scaling, PCA, or marker finding on large datasets, consider reducing the number of variable features (--n_variable_features), the number of PCA dimensions (--pca_dimensions), or processing samples individually rather than in parallel (or reduce the number of parallel jobs).

## Author: Vibhav Setlur
## Date: 2025-04-20

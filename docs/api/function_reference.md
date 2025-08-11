# Function Reference

This document provides a comprehensive reference for the functions and classes available in the snRNA-seq pipeline.

## Core Pipeline Functions

### `run_single_sample_pipeline()`

Main function for processing a single sample.

**Parameters:**
- `input_file` (string): Path to input file (H5 or RDS)
- `project_name` (string): Name for the analysis project
- `qc_params` (list): Quality control parameters
- `processing_params` (list): Processing parameters
- `clustering_params` (list): Clustering parameters
- `output_dir` (string): Output directory path
- `find_markers` (logical): Whether to find cluster markers
- `verbose` (logical): Enable verbose output

**Returns:**
- Seurat object with analysis results

**Example:**
```r
library(snrnaPipeline)

# Run single sample analysis
seurat_obj <- run_single_sample_pipeline(
  input_file = "data/sample.h5",
  project_name = "MyProject",
  qc_params = list(
    min_features = 200,
    max_mt_percent = 20,
    min_counts = 1000
  ),
  processing_params = list(
    n_variable_features = 2000,
    pca_dimensions = 15
  ),
  clustering_params = list(
    resolution = 0.5,
    algorithm = "leiden"
  ),
  find_markers = TRUE,
  verbose = TRUE
)
```

### `run_multi_sample_pipeline()`

Process multiple samples with parallel processing.

**Parameters:**
- `input_files` (vector): Paths to input files
- `project_name` (string): Name for the analysis project
- `sample_configs` (list): Sample-specific configurations
- `parallel` (logical): Enable parallel processing
- `n_cores` (integer): Number of cores to use
- `output_dir` (string): Output directory path

**Returns:**
- List of processed Seurat objects

**Example:**
```r
# Run multi-sample analysis
results <- run_multi_sample_pipeline(
  input_files = c("data/sample1.h5", "data/sample2.h5"),
  project_name = "MultiSampleProject",
  sample_configs = list(
    sample1 = list(
      qc = list(min_features = 200),
      clustering = list(resolution = 0.5)
    ),
    sample2 = list(
      qc = list(min_features = 300),
      clustering = list(resolution = 0.8)
    )
  ),
  parallel = TRUE,
  n_cores = 4
)
```

## Quality Control Functions

### `perform_quality_control()`

Perform quality control filtering on Seurat object.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `min_features` (integer): Minimum features per cell
- `max_features` (integer): Maximum features per cell
- `min_counts` (integer): Minimum UMI counts per cell
- `max_counts` (integer): Maximum UMI counts per cell
- `max_mt_percent` (numeric): Maximum mitochondrial percentage
- `min_cells` (integer): Minimum cells expressing a gene

**Returns:**
- Filtered Seurat object

**Example:**
```r
# Perform QC filtering
filtered_obj <- perform_quality_control(
  seurat_obj = raw_data,
  min_features = 200,
  max_features = 6000,
  min_counts = 1000,
  max_counts = 25000,
  max_mt_percent = 20,
  min_cells = 3
)
```

### `calculate_qc_metrics()`

Calculate quality control metrics for a Seurat object.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `species` (string): Species for mitochondrial gene detection

**Returns:**
- Data frame with QC metrics

**Example:**
```r
# Calculate QC metrics
qc_metrics <- calculate_qc_metrics(
  seurat_obj = raw_data,
  species = "human"
)

# View metrics
print(qc_metrics)
```

## Processing Functions

### `normalize_data()`

Normalize gene expression data.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `method` (string): Normalization method ("LogNormalize", "CLR", "RC")
- `scale_factor` (integer): Scale factor for normalization

**Returns:**
- Normalized Seurat object

**Example:**
```r
# Normalize data
normalized_obj <- normalize_data(
  seurat_obj = filtered_obj,
  method = "CLR",
  scale_factor = 10000
)
```

### `find_variable_features()`

Identify highly variable genes.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `n_features` (integer): Number of variable features to select
- `selection_method` (string): Feature selection method

**Returns:**
- Seurat object with variable features identified

**Example:**
```r
# Find variable features
var_features_obj <- find_variable_features(
  seurat_obj = normalized_obj,
  n_features = 2000,
  selection_method = "vst"
)
```

### `run_dimensionality_reduction()`

Perform PCA and UMAP dimensionality reduction.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `pca_dims` (integer): Number of PCA dimensions
- `umap_neighbors` (integer): Number of UMAP neighbors
- `umap_min_dist` (numeric): UMAP minimum distance

**Returns:**
- Seurat object with dimensionality reduction results

**Example:**
```r
# Run dimensionality reduction
reduced_obj <- run_dimensionality_reduction(
  seurat_obj = var_features_obj,
  pca_dims = 15,
  umap_neighbors = 30,
  umap_min_dist = 0.3
)
```

## Clustering Functions

### `perform_clustering()`

Perform graph-based clustering.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `resolution` (numeric): Clustering resolution
- `algorithm` (string): Clustering algorithm ("louvain", "leiden", "slm")
- `min_cluster_size` (integer): Minimum cluster size

**Returns:**
- Seurat object with clustering results

**Example:**
```r
# Perform clustering
clustered_obj <- perform_clustering(
  seurat_obj = reduced_obj,
  resolution = 0.5,
  algorithm = "leiden",
  min_cluster_size = 10
)
```

### `find_cluster_markers()`

Find marker genes for each cluster.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `min_pct` (numeric): Minimum percentage of cells expressing gene
- `logfc_threshold` (numeric): Log fold change threshold
- `test_use` (string): Statistical test to use

**Returns:**
- Data frame with marker genes

**Example:**
```r
# Find cluster markers
markers <- find_cluster_markers(
  seurat_obj = clustered_obj,
  min_pct = 0.25,
  logfc_threshold = 0.25,
  test_use = "wilcox"
)

# View top markers
head(markers)
```

## Visualization Functions

### `create_umap_plot()`

Create UMAP visualization.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `group_by` (string): Grouping variable
- `label` (logical): Whether to label clusters
- `colors` (vector): Color palette

**Returns:**
- ggplot object

**Example:**
```r
# Create UMAP plot
umap_plot <- create_umap_plot(
  seurat_obj = clustered_obj,
  group_by = "seurat_clusters",
  label = TRUE,
  colors = c("red", "blue", "green", "purple")
)

# Display plot
print(umap_plot)
```

### `create_feature_plot()`

Create feature plot for specific genes.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `features` (vector): Genes to plot
- `reduction` (string): Dimensionality reduction to use
- `ncol` (integer): Number of columns in plot grid

**Returns:**
- ggplot object

**Example:**
```r
# Create feature plot
feature_plot <- create_feature_plot(
  seurat_obj = clustered_obj,
  features = c("GENE1", "GENE2", "GENE3"),
  reduction = "umap",
  ncol = 3
)

# Display plot
print(feature_plot)
```

### `create_heatmap()`

Create heatmap of marker genes.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `markers` (data.frame): Marker genes data frame
- `top_n` (integer): Number of top markers per cluster
- `scale` (logical): Whether to scale the data

**Returns:**
- ComplexHeatmap object

**Example:**
```r
# Create heatmap
heatmap <- create_heatmap(
  seurat_obj = clustered_obj,
  markers = markers,
  top_n = 10,
  scale = TRUE
)

# Display heatmap
print(heatmap)
```

## Utility Functions

### `load_data()`

Load data from various file formats.

**Parameters:**
- `file_path` (string): Path to input file
- `file_type` (string): File type ("h5", "rds", "mtx")

**Returns:**
- Seurat object

**Example:**
```r
# Load H5 file
h5_data <- load_data("data/sample.h5", "h5")

# Load RDS file
rds_data <- load_data("data/sample.rds", "rds")
```

### `save_results()`

Save analysis results.

**Parameters:**
- `seurat_obj` (Seurat): Seurat object to save
- `output_dir` (string): Output directory
- `project_name` (string): Project name
- `save_plots` (logical): Whether to save plots

**Returns:**
- NULL

**Example:**
```r
# Save results
save_results(
  seurat_obj = clustered_obj,
  output_dir = "results",
  project_name = "MyProject",
  save_plots = TRUE
)
```

### `generate_report()`

Generate HTML report of analysis results.

**Parameters:**
- `seurat_obj` (Seurat): Seurat object
- `markers` (data.frame): Marker genes
- `output_file` (string): Output HTML file path

**Returns:**
- NULL

**Example:**
```r
# Generate report
generate_report(
  seurat_obj = clustered_obj,
  markers = markers,
  output_file = "results/analysis_report.html"
)
```

## Configuration Functions

### `load_config()`

Load configuration from YAML file.

**Parameters:**
- `config_file` (string): Path to configuration file

**Returns:**
- List with configuration parameters

**Example:**
```r
# Load configuration
config <- load_config("config/settings.yaml")

# Access parameters
qc_params <- config$quality_control
clustering_params <- config$clustering
```

### `validate_config()`

Validate configuration parameters.

**Parameters:**
- `config` (list): Configuration list

**Returns:**
- Logical indicating if configuration is valid

**Example:**
```r
# Validate configuration
is_valid <- validate_config(config)

if (!is_valid) {
  stop("Invalid configuration parameters")
}
```

## Multi-Sample Analysis Functions

### `integrate_samples()`

Integrate multiple samples using Seurat's integration methods.

**Parameters:**
- `sample_list` (list): List of Seurat objects
- `integration_method` (string): Integration method ("cca", "rpca", "harmony")
- `reference` (vector): Reference samples for integration

**Returns:**
- Integrated Seurat object

**Example:**
```r
# Integrate samples
integrated_obj <- integrate_samples(
  sample_list = list(sample1, sample2, sample3),
  integration_method = "cca",
  reference = c(1, 2)
)
```

### `compare_samples()`

Compare samples and generate comparison plots.

**Parameters:**
- `sample_list` (list): List of Seurat objects
- `comparison_type` (string): Type of comparison ("qc", "clustering", "markers")

**Returns:**
- List of comparison plots

**Example:**
```r
# Compare samples
comparison_plots <- compare_samples(
  sample_list = list(sample1, sample2),
  comparison_type = "qc"
)

# Display plots
for (plot in comparison_plots) {
  print(plot)
}
```

## Error Handling and Logging

### `setup_logging()`

Setup logging configuration.

**Parameters:**
- `log_file` (string): Log file path
- `log_level` (string): Logging level ("INFO", "WARNING", "ERROR")

**Returns:**
- NULL

**Example:**
```r
# Setup logging
setup_logging(
  log_file = "logs/pipeline.log",
  log_level = "INFO"
)
```

### `handle_errors()`

Error handling wrapper for pipeline functions.

**Parameters:**
- `expr` (expression): Expression to evaluate
- `error_msg` (string): Custom error message

**Returns:**
- Result of expression or error message

**Example:**
```r
# Handle errors
result <- handle_errors({
  run_single_sample_pipeline(input_file, project_name)
}, "Pipeline execution failed")
```

## Performance Functions

### `optimize_memory()`

Optimize memory usage for large datasets.

**Parameters:**
- `seurat_obj` (Seurat): Input Seurat object
- `memory_limit` (numeric): Memory limit in GB

**Returns:**
- Memory-optimized Seurat object

**Example:**
```r
# Optimize memory
optimized_obj <- optimize_memory(
  seurat_obj = large_dataset,
  memory_limit = 16
)
```

### `parallel_processing()`

Setup parallel processing for multi-core systems.

**Parameters:**
- `n_cores` (integer): Number of cores to use
- `memory_limit` (numeric): Memory limit per core

**Returns:**
- NULL

**Example:**
```r
# Setup parallel processing
parallel_processing(
  n_cores = 8,
  memory_limit = 4
)
```

This function reference provides comprehensive documentation for all major functions in the pipeline. For additional examples and use cases, refer to the vignettes and example workflows.

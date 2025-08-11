# Advanced Usage Guide

This guide covers advanced features and workflows for experienced users, including multi-sample analysis, custom parameter optimization, and integration with other computational biology tools.

## Multi-Sample Analysis

### Overview

Multi-sample analysis enables you to process multiple samples simultaneously with parallel processing and sample-specific parameter optimization. This is particularly useful for:

- Comparative studies across different conditions
- Batch processing of large datasets
- Integration of multiple experimental replicates
- Cross-sample quality control and normalization

### Basic Multi-Sample Workflow

```bash
# Process multiple H5 files
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample1.h5 data/sample2.h5 data/sample3.h5 \
  --project_name MultiSampleProject \
  --parallel \
  --n_cores 4 \
  --find_markers
```

### Sample-Specific Configuration

Create a YAML configuration file to customize parameters for individual samples:

```yaml
samples:
  control_sample:
    qc:
      min_features: 200
      max_mt_percent: 20
      min_counts: 1000
    processing:
      n_variable_features: 2000
      pca_dimensions: 15
    clustering:
      resolution: 0.5
      algorithm: "leiden"
  
  treatment_sample:
    qc:
      min_features: 300
      max_mt_percent: 15
      min_counts: 1500
    processing:
      n_variable_features: 2500
      pca_dimensions: 20
    clustering:
      resolution: 0.8
      algorithm: "louvain"
```

Run with custom configuration:

```bash
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/control.h5 data/treatment.h5 \
  --project_name TreatmentStudy \
  --sample_configs config/sample_configs.yaml \
  --parallel
```

### Output Structure

Multi-sample analysis creates an organized output structure:

```
MultiSampleProject_outputs/
├── individual_samples/
│   ├── sample1/
│   │   ├── sample1_processed.rds
│   │   ├── sample1_UMAP.png
│   │   ├── sample1_cluster_markers.csv
│   │   └── sample1_qc_summary.txt
│   ├── sample2/
│   │   └── ...
│   └── sample3/
│       └── ...
├── comparisons/
│   ├── MultiSampleProject_cell_counts_comparison.png
│   ├── MultiSampleProject_gene_counts_comparison.png
│   ├── MultiSampleProject_cluster_counts_comparison.png
│   └── MultiSampleProject_quality_metrics_comparison.png
├── combined_analysis/
│   ├── integrated_seurat_object.rds
│   ├── batch_corrected_umap.png
│   └── cross_sample_clustering.png
├── reports/
│   ├── MultiSampleProject_summary.html
│   └── MultiSampleProject_quality_report.html
└── MultiSampleProject_multi_sample_summary.csv
```

## Custom Workflows

### Parameter Optimization

#### Quality Control Optimization

```r
# Load your data
library(Seurat)
seurat_obj <- readRDS("data/sample.rds")

# Explore QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Determine optimal thresholds
qc_stats <- data.frame(
  n_features = seurat_obj$nFeature_RNA,
  n_counts = seurat_obj$nCount_RNA,
  mt_percent = seurat_obj$percent.mt
)

# Calculate percentiles for threshold selection
quantile(qc_stats$n_features, probs = c(0.01, 0.99))
quantile(qc_stats$n_counts, probs = c(0.01, 0.99))
quantile(qc_stats$mt_percent, probs = c(0.95, 0.99))
```

#### Clustering Resolution Optimization

```r
# Test multiple resolutions
resolutions <- c(0.1, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5)

# Run clustering at different resolutions
for(res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res, algorithm = 1)
  col_name <- paste0("RNA_snn_res.", res)
  print(paste("Resolution", res, ":", length(unique(seurat_obj@meta.data[[col_name]]))))
}

# Plot UMAP for different resolutions
for(res in resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = col_name, label = TRUE)
  print(p + ggtitle(paste("Resolution", res)))
}
```

### Custom Analysis Scripts

#### Integration with External Tools

```r
# Example: Integration with Monocle3 for trajectory analysis
library(monocle3)

# Convert Seurat object to Monocle3
cds <- as.cell_data_set(seurat_obj)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

# Plot trajectory
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = FALSE)
```

#### Custom Visualization

```r
# Create publication-ready plots
library(ggplot2)
library(patchwork)

# Custom UMAP with specific colors
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  scale_color_manual(values = c("red", "blue", "green", "purple")) +
  theme_minimal() +
  ggtitle("Cell Clusters")

# Feature plot with custom settings
feature_plot <- FeaturePlot(seurat_obj, features = c("GENE1", "GENE2"), 
                           reduction = "umap", ncol = 2) &
  theme_minimal() &
  theme(legend.position = "bottom")

# Combine plots
combined_plot <- umap_plot / feature_plot
ggsave("custom_visualization.png", combined_plot, width = 12, height = 8)
```

## Integration with Other Tools

### Cell Ranger Integration

```bash
# Process Cell Ranger output directly
Rscript scripts/run_pipeline_terminal.R \
  --h5_input /path/to/cellranger/outs/filtered_feature_bc_matrix.h5 \
  --project_name CellRangerAnalysis
```

### SouporCell Doublet Detection

```bash
# Run with SouporCell doublet detection
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample1.h5 data/sample2.h5 \
  --soupor_cell_doublet_inputs data/soupor1.tsv data/soupor2.tsv \
  --project_name DoubletFilteredAnalysis
```

### Integration with Seurat Workflows

```r
# Load pipeline output into existing Seurat workflow
pipeline_output <- readRDS("MyProject_outputs/individual_samples/sample/sample_processed.rds")

# Continue with custom analysis
pipeline_output <- FindClusters(pipeline_output, resolution = 0.8)
pipeline_output <- RunUMAP(pipeline_output, dims = 1:30)

# Add custom metadata
pipeline_output$custom_annotation <- "your_annotation"

# Save updated object
saveRDS(pipeline_output, "updated_analysis.rds")
```

## Performance Optimization

### Memory Management

```bash
# Increase memory allocation
export R_MAX_MEM_SIZE=32G

# Use parallel processing
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/*.h5 \
  --project_name LargeDataset \
  --parallel \
  --n_cores 8
```

### Batch Processing

Create a batch processing script:

```bash
#!/bin/bash
# batch_process.sh

for file in data/*.h5; do
    sample_name=$(basename "$file" .h5)
    echo "Processing $sample_name..."
    
    Rscript scripts/run_pipeline_terminal.R \
      --h5_input "$file" \
      --project_name "$sample_name" \
      --find_markers \
      --verbose
done
```

## Quality Control Best Practices

### Data Quality Assessment

```r
# Comprehensive QC assessment
qc_assessment <- function(seurat_obj) {
  # Calculate QC metrics
  qc_metrics <- data.frame(
    total_cells = ncol(seurat_obj),
    total_genes = nrow(seurat_obj),
    median_features = median(seurat_obj$nFeature_RNA),
    median_counts = median(seurat_obj$nCount_RNA),
    median_mt_percent = median(seurat_obj$percent.mt),
    cells_after_qc = sum(seurat_obj$nFeature_RNA >= 200 & 
                        seurat_obj$nCount_RNA >= 1000 & 
                        seurat_obj$percent.mt <= 20)
  )
  
  return(qc_metrics)
}
```

### Cross-Sample Quality Control

```r
# Compare quality across samples
compare_quality <- function(sample_list) {
  quality_comparison <- data.frame()
  
  for(sample in names(sample_list)) {
    qc <- qc_assessment(sample_list[[sample]])
    qc$sample <- sample
    quality_comparison <- rbind(quality_comparison, qc)
  }
  
  return(quality_comparison)
}
```

## Troubleshooting Advanced Issues

### Memory Issues

```bash
# Monitor memory usage
htop

# Use swap space for large datasets
sudo fallocate -l 32G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

### Parallel Processing Issues

```r
# Check available cores
library(parallel)
detectCores()

# Set number of cores manually
options(mc.cores = 4)
```

### Integration Issues

```r
# Check package compatibility
sessionInfo()

# Update packages if needed
update.packages(ask = FALSE)
```

## Example Advanced Workflow

Here's a complete advanced workflow for a multi-condition study:

```bash
#!/bin/bash
# advanced_workflow.sh

# 1. Setup
export R_MAX_MEM_SIZE=32G

# 2. Process control samples
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/control_*.h5 \
  --project_name ControlSamples \
  --sample_configs config/control_config.yaml \
  --parallel \
  --n_cores 4 \
  --find_markers

# 3. Process treatment samples
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/treatment_*.h5 \
  --project_name TreatmentSamples \
  --sample_configs config/treatment_config.yaml \
  --parallel \
  --n_cores 4 \
  --find_markers

# 4. Combined analysis
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/*.h5 \
  --project_name CombinedAnalysis \
  --sample_configs config/combined_config.yaml \
  --parallel \
  --n_cores 8 \
  --find_markers
```

This workflow provides a comprehensive analysis suitable for publication and further downstream analysis.

# Example Workflows

This document provides practical examples of how to use the snRNA-seq pipeline for different research scenarios.

## Basic Single Sample Analysis

### Example 1: Brain Tissue Analysis

**Scenario**: Analyze single-nucleus RNA-seq data from brain tissue to identify cell types.

```bash
# Setup
git clone <repository-url>
cd snRNA-seq-Pipeline
./launch.sh setup

# Run analysis
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/brain_sample.h5 \
  --project_name BrainAnalysis \
  --min_features 200 \
  --max_mt_percent 20 \
  --n_variable_features 2000 \
  --clustering_resolution 0.5 \
  --find_markers \
  --verbose
```

**Expected Output**:
- Cell type clusters with marker genes
- UMAP visualization of cell types
- Quality control metrics
- Differential expression analysis

### Example 2: Cancer Tissue Analysis

**Scenario**: Analyze tumor tissue with more stringent quality control.

```bash
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/tumor_sample.h5 \
  --project_name TumorAnalysis \
  --min_features 300 \
  --min_counts 1500 \
  --max_mt_percent 15 \
  --n_variable_features 2500 \
  --clustering_resolution 0.8 \
  --clustering_algorithm louvain \
  --find_markers
```

## Multi-Sample Comparative Analysis

### Example 3: Control vs Treatment Study

**Scenario**: Compare gene expression between control and treatment conditions.

**Step 1: Create sample configuration**
```yaml
# config/treatment_study.yaml
samples:
  control:
    qc:
      min_features: 200
      max_mt_percent: 20
    clustering:
      resolution: 0.5
      algorithm: "leiden"
  
  treatment:
    qc:
      min_features: 200
      max_mt_percent: 20
    clustering:
      resolution: 0.5
      algorithm: "leiden"
```

**Step 2: Run analysis**
```bash
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/control.h5 data/treatment.h5 \
  --project_name TreatmentStudy \
  --sample_configs config/treatment_study.yaml \
  --parallel \
  --n_cores 4 \
  --find_markers
```

**Step 3: Analyze results**
```r
# Load results
library(Seurat)

# Load individual samples
control <- readRDS("TreatmentStudy_outputs/individual_samples/control/control_processed.rds")
treatment <- readRDS("TreatmentStudy_outputs/individual_samples/treatment/treatment_processed.rds")

# Load integrated analysis
integrated <- readRDS("TreatmentStudy_outputs/combined_analysis/integrated_seurat_object.rds")

# Compare cell type proportions
control_props <- table(control$seurat_clusters) / ncol(control)
treatment_props <- table(treatment$seurat_clusters) / ncol(treatment)

# Statistical test for differences
chisq.test(rbind(control_props, treatment_props))
```

### Example 4: Time Course Analysis

**Scenario**: Analyze samples collected at different time points.

```bash
# Process all time points
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/day0.h5 data/day3.h5 data/day7.h5 data/day14.h5 \
  --project_name TimeCourse \
  --sample_configs config/timecourse.yaml \
  --parallel \
  --n_cores 4
```

**Configuration for time course**:
```yaml
# config/timecourse.yaml
samples:
  day0:
    qc:
      min_features: 200
    clustering:
      resolution: 0.6
  day3:
    qc:
      min_features: 200
    clustering:
      resolution: 0.6
  day7:
    qc:
      min_features: 200
    clustering:
      resolution: 0.6
  day14:
    qc:
      min_features: 200
    clustering:
      resolution: 0.6
```

## Large Dataset Analysis

### Example 5: High-Throughput Screening

**Scenario**: Process multiple samples from a high-throughput experiment.

```bash
#!/bin/bash
# batch_process.sh

# Process all samples in parallel
for file in data/screen_*.h5; do
    sample_name=$(basename "$file" .h5)
    echo "Processing $sample_name..."
    
    Rscript scripts/run_pipeline_terminal.R \
      --h5_input "$file" \
      --project_name "$sample_name" \
      --min_features 200 \
      --n_variable_features 2000 \
      --clustering_resolution 0.5 \
      --find_markers \
      --output_format png &
done

# Wait for all processes to complete
wait

echo "All samples processed!"
```

### Example 6: Memory-Optimized Analysis

**Scenario**: Analyze large datasets with limited memory.

```bash
# Use conservative settings for large datasets
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/large_sample.h5 \
  --project_name LargeDataset \
  --min_features 500 \
  --min_counts 2000 \
  --n_variable_features 1500 \
  --pca_dimensions 20 \
  --save_intermediate FALSE \
  --output_format png
```

## Integration with External Tools

### Example 7: Cell Ranger Integration

**Scenario**: Process Cell Ranger output directly.

```bash
# Process Cell Ranger output
Rscript scripts/run_pipeline_terminal.R \
  --h5_input /path/to/cellranger/outs/filtered_feature_bc_matrix.h5 \
  --project_name CellRangerAnalysis \
  --min_features 200 \
  --max_mt_percent 20
```

### Example 8: SouporCell Doublet Detection

**Scenario**: Remove doublets using SouporCell results.

```bash
# Run with doublet detection
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample1.h5 data/sample2.h5 \
  --soupor_cell_doublet_inputs data/soupor1.tsv data/soupor2.tsv \
  --project_name DoubletFiltered \
  --parallel
```

## Custom Analysis Workflows

### Example 9: Trajectory Analysis

**Scenario**: Perform trajectory analysis after pipeline processing.

```r
# Load pipeline output
library(Seurat)
library(monocle3)

# Load processed data
seurat_obj <- readRDS("MyProject_outputs/individual_samples/sample/sample_processed.rds")

# Convert to Monocle3
cds <- as.cell_data_set(seurat_obj)

# Preprocess
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

# Plot trajectory
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = FALSE)
```

### Example 10: Custom Visualization

**Scenario**: Create publication-ready figures.

```r
# Load data
library(Seurat)
library(ggplot2)
library(patchwork)

seurat_obj <- readRDS("MyProject_outputs/individual_samples/sample/sample_processed.rds")

# Create custom UMAP
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) +
  theme_minimal() +
  ggtitle("Cell Type Clusters") +
  theme(legend.position = "bottom")

# Create feature plots
feature_plot <- FeaturePlot(seurat_obj, 
                           features = c("GENE1", "GENE2", "GENE3"), 
                           reduction = "umap", 
                           ncol = 3) &
  theme_minimal() &
  theme(legend.position = "bottom")

# Combine plots
combined_plot <- umap_plot / feature_plot

# Save publication-ready figure
ggsave("publication_figure.png", combined_plot, 
       width = 12, height = 8, dpi = 300)
```

## Quality Control Workflows

### Example 11: QC Parameter Optimization

**Scenario**: Determine optimal QC parameters for your data.

```r
# Load raw data
library(Seurat)
seurat_obj <- readRDS("data/raw_sample.rds")

# Analyze QC distributions
qc_stats <- data.frame(
  n_features = seurat_obj$nFeature_RNA,
  n_counts = seurat_obj$nCount_RNA,
  mt_percent = seurat_obj$percent.mt
)

# Calculate percentiles
percentiles <- data.frame(
  metric = c("n_features", "n_counts", "mt_percent"),
  p1 = c(quantile(qc_stats$n_features, 0.01),
          quantile(qc_stats$n_counts, 0.01),
          quantile(qc_stats$mt_percent, 0.95)),
  p5 = c(quantile(qc_stats$n_features, 0.05),
          quantile(qc_stats$n_counts, 0.05),
          quantile(qc_stats$mt_percent, 0.99)),
  p95 = c(quantile(qc_stats$n_features, 0.95),
           quantile(qc_stats$n_counts, 0.95),
           quantile(qc_stats$mt_percent, 0.01)),
  p99 = c(quantile(qc_stats$n_features, 0.99),
           quantile(qc_stats$n_counts, 0.99),
           quantile(qc_stats$mt_percent, 0.05))
)

print(percentiles)
```

### Example 12: Cross-Sample QC

**Scenario**: Compare quality metrics across multiple samples.

```r
# Load multiple samples
samples <- list()
sample_files <- list.files("MultiProject_outputs/individual_samples/", 
                          full.names = TRUE)

for(sample_dir in sample_files) {
  sample_name <- basename(sample_dir)
  rds_file <- file.path(sample_dir, paste0(sample_name, "_processed.rds"))
  if(file.exists(rds_file)) {
    samples[[sample_name]] <- readRDS(rds_file)
  }
}

# Compare QC metrics
qc_comparison <- data.frame()
for(sample_name in names(samples)) {
  seurat_obj <- samples[[sample_name]]
  qc_metrics <- data.frame(
    sample = sample_name,
    total_cells = ncol(seurat_obj),
    median_features = median(seurat_obj$nFeature_RNA),
    median_counts = median(seurat_obj$nCount_RNA),
    median_mt_percent = median(seurat_obj$percent.mt)
  )
  qc_comparison <- rbind(qc_comparison, qc_metrics)
}

print(qc_comparison)
```

## Performance Optimization Examples

### Example 13: Parallel Processing Setup

**Scenario**: Optimize processing for multi-core systems.

```bash
#!/bin/bash
# parallel_processing.sh

# Set environment variables
export R_MAX_MEM_SIZE=32G
export OMP_NUM_THREADS=8

# Process with maximum parallelization
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample*.h5 \
  --project_name ParallelAnalysis \
  --parallel \
  --n_cores 8 \
  --find_markers
```

### Example 14: Memory-Efficient Processing

**Scenario**: Process large datasets with limited memory.

```bash
# Use conservative memory settings
export R_MAX_MEM_SIZE=8G

# Process with memory-efficient settings
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/large_sample.h5 \
  --project_name MemoryEfficient \
  --min_features 500 \
  --n_variable_features 1500 \
  --pca_dimensions 15 \
  --save_intermediate FALSE
```

These examples provide practical starting points for different research scenarios. Modify parameters and workflows based on your specific data characteristics and research questions.

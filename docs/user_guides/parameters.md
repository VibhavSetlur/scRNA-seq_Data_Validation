# Parameter Reference

This document provides a complete reference for all parameters available in the snRNA-seq pipeline.

## Command Line Parameters

### Input Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--h5_input` | string | - | Path to 10X Genomics H5 file |
| `--rds_input` | string | - | Path to Seurat RDS file |
| `--project_name` | string | "SeuratProject" | Name for the analysis project |
| `--working_dir` | string | "." | Working directory for analysis |
| `--soupor_cell_doublet_input` | string | - | Path to SouporCell doublet detection output |

### Quality Control Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--min_features` | integer | 200 | Minimum features per cell |
| `--max_features` | integer | 6000 | Maximum features per cell |
| `--min_counts` | integer | 1000 | Minimum UMI counts per cell |
| `--max_counts` | integer | 25000 | Maximum UMI counts per cell |
| `--max_mt_percent` | float | 20.0 | Maximum mitochondrial percentage |
| `--min_cells` | integer | 3 | Minimum cells expressing a gene |

### Processing Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--n_variable_features` | integer | 2000 | Number of highly variable genes |
| `--normalization_method` | string | "CLR" | Normalization method (LogNormalize\|RC\|CLR) |
| `--scaling_method` | string | "negbinom" | Scaling method (negbinom\|linear) |
| `--pca_dimensions` | integer | 15 | Number of PCA dimensions |
| `--scale_factor` | integer | 10000 | Scale factor for normalization |

### Clustering Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--clustering_resolution` | float | 0.5 | Clustering resolution |
| `--clustering_algorithm` | string | "leiden" | Clustering algorithm (louvain\|multilevel\|leiden\|slm) |
| `--min_cluster_size` | integer | 10 | Minimum cluster size |
| `--umap_n_neighbors` | integer | 30 | Number of UMAP neighbors |
| `--umap_min_dist` | float | 0.3 | UMAP minimum distance |

### Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--find_markers` | flag | FALSE | Perform differential expression analysis |
| `--save_intermediate` | flag | FALSE | Save intermediate results |
| `--random_seed` | integer | 42 | Random seed for reproducibility |
| `--verbose` | flag | FALSE | Enable verbose output |
| `--output_format` | string | "png" | Output format for plots (png\|pdf\|svg) |

## Multi-Sample Parameters

### Input Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--h5_inputs` | string[] | - | Paths to multiple H5 files |
| `--rds_inputs` | string[] | - | Paths to multiple RDS files |
| `--soupor_cell_doublet_inputs` | string[] | - | Paths to multiple SouporCell outputs |
| `--sample_configs` | string | - | Path to sample-specific configuration file |

### Parallel Processing Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--parallel` | flag | TRUE | Enable parallel processing |
| `--n_cores` | integer | auto | Number of cores to use |

## Sample Configuration File

For multi-sample analysis, you can create a YAML configuration file to customize parameters for individual samples:

```yaml
samples:
  sample1:
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
      min_cluster_size: 15
  
  sample2:
    qc:
      min_features: 200
      max_mt_percent: 25
    processing:
      n_variable_features: 3000
    clustering:
      resolution: 0.3
      algorithm: "leiden"
```

### Configuration File Parameters

#### Quality Control Section
- `min_features`: Minimum features per cell
- `max_features`: Maximum features per cell
- `min_counts`: Minimum UMI counts per cell
- `max_counts`: Maximum UMI counts per cell
- `max_mt_percent`: Maximum mitochondrial percentage
- `min_cells`: Minimum cells expressing a gene

#### Processing Section
- `n_variable_features`: Number of highly variable genes
- `normalization_method`: Normalization method
- `scaling_method`: Scaling method
- `pca_dimensions`: Number of PCA dimensions
- `scale_factor`: Scale factor for normalization

#### Clustering Section
- `resolution`: Clustering resolution
- `algorithm`: Clustering algorithm
- `min_cluster_size`: Minimum cluster size
- `umap_n_neighbors`: Number of UMAP neighbors
- `umap_min_dist`: UMAP minimum distance

## Parameter Recommendations

### Quality Control Thresholds

**For high-quality data:**
- `min_features`: 200-300
- `max_features`: 6000-8000
- `min_counts`: 1000-2000
- `max_counts`: 25000-50000
- `max_mt_percent`: 15-25

**For lower quality data:**
- `min_features`: 100-200
- `max_features`: 8000-10000
- `min_counts`: 500-1000
- `max_counts`: 50000-100000
- `max_mt_percent`: 25-35

### Clustering Parameters

**For fine-grained clustering:**
- `resolution`: 0.8-1.2
- `algorithm`: "leiden" or "louvain"

**For broad clustering:**
- `resolution`: 0.3-0.5
- `algorithm`: "louvain" or "slm"

### Processing Parameters

**For standard analysis:**
- `n_variable_features`: 2000-3000
- `pca_dimensions`: 15-30
- `normalization_method`: "CLR"

**For specialized analysis:**
- `n_variable_features`: 1000-5000 (depending on data complexity)
- `pca_dimensions`: 10-50 (depending on data size)
- `normalization_method`: "LogNormalize" (for standard workflows)

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `SHINY_HOST` | "0.0.0.0" | Host address for Shiny app |
| `SHINY_PORT` | "3838" | Port for Shiny app |
| `R_MAX_MEM_SIZE` | "8G" | Maximum R memory allocation |

## Examples

### Basic Single Sample Analysis

```bash
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name MyProject \
  --min_features 200 \
  --max_mt_percent 20 \
  --n_variable_features 2000 \
  --clustering_resolution 0.5 \
  --find_markers
```

### Multi-Sample Analysis with Custom Config

```bash
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample1.h5 data/sample2.h5 \
  --project_name MultiProject \
  --sample_configs config/sample_configs.yaml \
  --parallel \
  --n_cores 4 \
  --find_markers
```

### High-Resolution Clustering

```bash
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name HighResProject \
  --clustering_resolution 1.0 \
  --clustering_algorithm leiden \
  --min_cluster_size 5 \
  --n_variable_features 3000 \
  --pca_dimensions 25
```

### Conservative Quality Control

```bash
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name ConservativeQC \
  --min_features 300 \
  --min_counts 1500 \
  --max_mt_percent 15 \
  --max_features 5000 \
  --max_counts 20000
```

# Parameter Reference

This document provides detailed descriptions of all parameters available in the snRNA-seq Pipeline.

## Data Input Parameters

### Input Type
- **H5 File**: 10X Genomics filtered feature-barcode matrix
- **RDS File**: Pre-processed Seurat object

### Project Settings
- **Project Name**: Name for your analysis project
- **Working Directory**: Directory for output files
- **SouporCell File**: Optional doublet detection results

## Quality Control Parameters

### Cell Filtering
- **Minimum Features**: Minimum number of genes detected per cell
  - Range: 0-10000
  - Default: 200
  - Typical: 200-500

- **Minimum Counts**: Minimum number of UMIs per cell
  - Range: 0-100000
  - Default: 1000
  - Typical: 1000-2000

- **Maximum Features**: Maximum number of genes detected per cell
  - Range: 0-50000
  - Default: 6000
  - Typical: 5000-8000

- **Maximum Counts**: Maximum number of UMIs per cell
  - Range: 0-1000000
  - Default: 25000
  - Typical: 20000-50000

- **Maximum Mitochondrial %**: Maximum percentage of mitochondrial genes
  - Range: 0-100
  - Default: 20
  - Typical: 10-25 (depends on tissue type)

### Gene Filtering
- **Minimum Cells**: Minimum number of cells expressing a gene
  - Range: 0-100
  - Default: 3
  - Typical: 3-10

## Processing Parameters

### Normalization
- **Normalization Method**:
  - **CLR**: Centered Log Ratio (recommended for UMI data)
  - **LogNormalize**: Log normalization
  - **RC**: Relative Counts

- **Scale Factor**: Factor for normalization
  - Range: 1000-100000
  - Default: 10000
  - Typical: 10000

### Feature Selection
- **Number of Variable Features**: Number of highly variable genes
  - Range: 100-10000
  - Default: 2000
  - Typical: 2000-3000

### Scaling
- **Scaling Method**:
  - **negbinom**: Negative binomial (recommended for UMI data)
  - **linear**: Linear scaling

### Dimensionality Reduction
- **PCA Dimensions**: Number of principal components
  - Range: 1-100
  - Default: 15
  - Typical: 15-30

## Clustering Parameters

### Clustering Algorithm
- **leiden**: Leiden algorithm (recommended)
- **louvain**: Louvain algorithm
- **multilevel**: Multilevel algorithm
- **slm**: SLM algorithm

### Clustering Settings
- **Clustering Resolution**: Controls number of clusters
  - Range: 0.1-2.0
  - Default: 0.5
  - Typical: 0.3-1.0

- **Minimum Cluster Size**: Minimum cells per cluster
  - Range: 1-1000
  - Default: 10
  - Typical: 10-50

### UMAP Parameters
- **UMAP n_neighbors**: Number of neighbors for UMAP
  - Range: 5-100
  - Default: 30
  - Typical: 30-50

- **UMAP min_dist**: Minimum distance between points
  - Range: 0.01-1.0
  - Default: 0.3
  - Typical: 0.1-0.5

## Visualization Parameters

### Plot Settings
- **Plot Theme**: Visual theme for plots
  - **classic**: Classic ggplot theme
  - **minimal**: Minimal theme
  - **bw**: Black and white theme

- **Color Palette**: Color scheme for plots
  - **viridis**: Viridis color palette (recommended)
  - **rainbow**: Rainbow color palette
  - **brewer**: ColorBrewer palettes

### Plot Elements
- **Point Size**: Size of points in scatter plots
  - Range: 0.1-5.0
  - Default: 0.7
  - Typical: 0.5-1.0

- **Label Size**: Size of text labels
  - Range: 1-10
  - Default: 3
  - Typical: 2-5

- **Title Size**: Size of plot titles
  - Range: 8-24
  - Default: 14
  - Typical: 12-16

### Output Format
- **Output Format**: File format for plots
  - **png**: Portable Network Graphics
  - **pdf**: Portable Document Format
  - **svg**: Scalable Vector Graphics

## Analysis Options

### Marker Analysis
- **Find Cluster Markers**: Perform differential expression analysis
  - Default: TRUE
  - Generates marker gene lists for each cluster

### Output Options
- **Save Intermediate Results**: Save intermediate processing steps
  - Default: TRUE
  - Useful for debugging and step-by-step analysis

- **Verbose Output**: Show detailed progress information
  - Default: TRUE
  - Helpful for monitoring long analyses

### Reproducibility
- **Random Seed**: Seed for random number generation
  - Range: 1-10000
  - Default: 42
  - Ensures reproducible results

## Parameter Recommendations

### For High-Quality Data
- Lower QC thresholds
- Higher variable features
- More PCA dimensions

### For Low-Quality Data
- Higher QC thresholds
- Fewer variable features
- Fewer PCA dimensions

### For Large Datasets
- Higher minimum features/counts
- More PCA dimensions
- Higher clustering resolution

### For Small Datasets
- Lower minimum features/counts
- Fewer PCA dimensions
- Lower clustering resolution

# Getting Started Guide

This guide will walk you through setting up and running your first single-nucleus RNA-seq analysis using the pipeline.

## Prerequisites

- **R** (version 4.0 or higher)
- **RStudio** (recommended for interactive use)
- **Git**
- **Docker** (optional, for containerized deployment)

## Installation

### Step 1: Clone the Repository

```bash
git clone <repository-url>
cd snRNA-seq-Pipeline
```

### Step 2: Setup Environment

```bash
# Run the automated setup script
./launch.sh setup
```

This script will:
- Install required R packages
- Set up configuration files
- Create necessary directories
- Verify the installation

### Step 3: Verify Installation

```bash
# Test the installation
Rscript scripts/run_pipeline_terminal.R --help
```

## Your First Analysis

### Option 1: Web Interface (Recommended for Beginners)

1. **Launch the Shiny app**:
   ```bash
   ./launch.sh shiny
   ```

2. **Open your browser** and navigate to `http://localhost:3838`

3. **Upload your data**:
   - Click "Browse" to select your H5 or RDS file
   - Supported formats: 10X Genomics H5 files, Seurat RDS objects

4. **Configure parameters**:
   - **Project Name**: Enter a descriptive name for your analysis
   - **Quality Control**: Adjust filtering thresholds if needed
   - **Processing**: Set number of variable features and PCA dimensions
   - **Clustering**: Choose clustering algorithm and resolution

5. **Run analysis**:
   - Click "Run Pipeline" to start the analysis
   - Monitor progress in the status panel

6. **Explore results**:
   - View QC plots and statistics
   - Examine clustering results
   - Download processed data and plots

### Option 2: Command Line Interface

1. **Prepare your data**:
   ```bash
   # Create a data directory
   mkdir -p data
   
   # Copy your data files
   cp /path/to/your/sample.h5 data/
   ```

2. **Run basic analysis**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --project_name MyFirstAnalysis
   ```

3. **Run with custom parameters**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --project_name MyFirstAnalysis \
     --min_features 300 \
     --max_mt_percent 15 \
     --clustering_resolution 0.8 \
     --find_markers
   ```

## Understanding Your Results

### Output Structure

After running the pipeline, you'll find your results in:

```
MyFirstAnalysis_outputs/
├── individual_samples/
│   └── sample/
│       ├── sample_processed.rds          # Processed Seurat object
│       ├── sample_UMAP.png               # UMAP visualization
│       ├── sample_cluster_markers.csv    # Cluster markers
│       └── sample_qc_summary.txt         # QC statistics
├── reports/
│   └── MyFirstAnalysis_summary.html      # Analysis report
└── MyFirstAnalysis_session_info.txt      # Analysis metadata
```

### Key Output Files

- **`*_processed.rds`**: Processed Seurat object with all analysis results
- **`*_UMAP.png`**: UMAP plot showing cell clusters
- **`*_cluster_markers.csv`**: Differential expression results for each cluster
- **`*_qc_summary.txt`**: Quality control statistics and filtering summary

### Loading Results in R

```r
# Load the processed Seurat object
seurat_obj <- readRDS("MyFirstAnalysis_outputs/individual_samples/sample/sample_processed.rds")

# View available metadata
colnames(seurat_obj@meta.data)

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# View cluster markers
markers <- read.csv("MyFirstAnalysis_outputs/individual_samples/sample/sample_cluster_markers.csv")
head(markers)
```

## Next Steps

Once you're comfortable with basic analysis:

1. **Multi-sample analysis**: Process multiple samples simultaneously
2. **Custom workflows**: Modify parameters for your specific research questions
3. **Advanced visualization**: Create publication-ready figures
4. **Integration**: Combine with other single-cell analysis tools

## Common Issues

### Installation Problems

**R package installation fails**:
```bash
# Update R packages
Rscript -e "update.packages(ask = FALSE)"

# Reinstall dependencies
Rscript setup/install_dependencies.R
```

**Permission errors**:
```bash
# Fix file permissions
chmod +x launch.sh
chmod +x scripts/*.R
```

### Runtime Issues

**Memory errors**:
```bash
# Increase R memory limit
export R_MAX_MEM_SIZE=16G
./launch.sh shiny
```

**Port conflicts**:
```bash
# Use different port
SHINY_PORT=8080 ./launch.sh shiny
```

## Getting Help

- Check the [troubleshooting guide](troubleshooting.md) for common issues
- Review [example workflows](../vignettes/) for advanced usage
- Open an issue on the repository for bugs or feature requests

## Example Workflow

Here's a complete example workflow for analyzing a 10X Genomics dataset:

```bash
# 1. Setup
git clone <repository-url>
cd snRNA-seq-Pipeline
./launch.sh setup

# 2. Prepare data
mkdir -p data
cp /path/to/10x_data.h5 data/

# 3. Run analysis
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/10x_data.h5 \
  --project_name BrainTissue \
  --min_features 200 \
  --max_mt_percent 20 \
  --n_variable_features 2000 \
  --clustering_resolution 0.5 \
  --find_markers \
  --verbose

# 4. Explore results
ls -la BrainTissue_outputs/
```

This workflow will generate a complete analysis with quality control, clustering, and marker identification suitable for publication.

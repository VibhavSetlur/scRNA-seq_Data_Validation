# snRNA-seq Analysis Pipeline

A comprehensive single-nucleus RNA sequencing analysis pipeline designed for computational biology research. Built with Seurat and R, this tool provides both interactive web-based analysis and command-line processing for flexible research workflows.

## Overview

This pipeline enables researchers to perform end-to-end single-nucleus RNA-seq analysis, from raw data processing to advanced clustering and marker identification. It supports both individual sample analysis and multi-sample comparative studies with parallel processing capabilities.

## Key Features

- **Multi-Sample Analysis**: Process multiple samples simultaneously with sample-specific parameter optimization
- **Flexible Input Formats**: Support for 10X Genomics H5 files and Seurat RDS objects
- **Quality Control**: Comprehensive filtering with customizable thresholds for mitochondrial content, feature counts, and cell quality
- **Advanced Clustering**: Multiple clustering algorithms (Louvain, Leiden, SLM) with resolution optimization
- **Marker Identification**: Automated differential expression analysis for cluster characterization
- **Interactive Visualization**: Web-based Shiny interface for exploratory data analysis
- **Batch Processing**: Command-line interface for high-throughput analysis
- **Reproducible Workflows**: Containerized deployment with version-controlled environments

## Quick Start

### Installation

```bash
# Clone repository
git clone <repository-url>
cd snRNA-seq-Pipeline

# Setup environment
./launch.sh setup
```

### Basic Usage

**Web Interface (Recommended for exploration):**
```bash
./launch.sh shiny
```

**Command Line (Recommended for batch processing):**
```bash
# Single sample analysis
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name MyProject

# Multi-sample analysis
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample1.h5 data/sample2.h5 \
  --project_name MultiSampleProject \
  --parallel
```

## Input Data Requirements

### Supported Formats
- **10X Genomics H5 files**: Raw count matrices from Cell Ranger output
- **Seurat RDS objects**: Pre-processed Seurat objects
- **SouporCell output**: Optional doublet detection results

### Data Quality Recommendations
- Minimum 200 features per cell
- Maximum 20% mitochondrial content
- Minimum 1,000 UMIs per cell
- Sample size: 1,000-100,000 cells per sample

## Analysis Workflow

1. **Data Loading**: Import H5 or RDS files with automatic format detection
2. **Quality Control**: Filter cells based on feature counts, UMI counts, and mitochondrial percentage
3. **Normalization**: Log-normalization or centered log-ratio transformation
4. **Feature Selection**: Identify highly variable genes
5. **Dimensionality Reduction**: PCA followed by UMAP/tSNE
6. **Clustering**: Graph-based clustering with multiple algorithm options
7. **Marker Analysis**: Differential expression analysis for cluster characterization
8. **Visualization**: Generate publication-ready plots and summary statistics

## Output Structure

```
ProjectName_outputs/
├── individual_samples/          # Per-sample results
├── combined_analysis/           # Integrated analysis
├── comparisons/                 # Cross-sample comparisons
├── reports/                     # Summary reports
└── ProjectName_session_info.txt # Analysis metadata
```

## Documentation

- **[User Guide](docs/user_guides/getting_started.md)**: Step-by-step tutorial for first-time users
- **[Parameter Reference](docs/user_guides/parameters.md)**: Complete parameter documentation
- **[Advanced Usage](docs/user_guides/advanced_usage.md)**: Multi-sample analysis and custom workflows
- **[Troubleshooting](docs/user_guides/troubleshooting.md)**: Common issues and solutions
- **[API Documentation](docs/api/)**: Function reference and examples

## Citation

If you use this pipeline in your research, please cite:

```
snRNA-seq Analysis Pipeline v1.0
Single-nucleus RNA sequencing analysis pipeline
https://github.com/your-repo/snRNA-seq-Pipeline
```

## Support

For questions and issues:
- Check the [troubleshooting guide](docs/user_guides/troubleshooting.md)
- Review [example workflows](docs/vignettes/)
- Open an issue on the repository

## License

MIT License - see [LICENSE](LICENSE) file for details.

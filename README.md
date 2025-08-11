# snRNA-seq Pipeline

A comprehensive single-nucleus RNA sequencing analysis pipeline built with Seurat, featuring both a web-based Shiny interface and command-line tools for flexible deployment and execution.

## Features

- **Multi-Sample Analysis**: Process multiple samples simultaneously with parallel processing
- **Sample-Specific Configurations**: Customize parameters for individual samples
- **Parallel Processing**: Efficient multi-core processing for faster analysis
- **Web Interface**: User-friendly Shiny application for interactive analysis
- **Command Line**: Terminal interface for batch processing and server deployment
- **Containerized**: Docker support for easy deployment and scaling
- **Scalable**: Designed for both local development and production server deployment
- **Flexible Input**: Supports both H5 files (10X Genomics) and RDS files (Seurat objects)
- **Quality Control**: Comprehensive QC filtering with customizable parameters
- **Advanced Analysis**: Clustering, dimensionality reduction, and marker identification
- **Visualization**: Automated generation of publication-ready plots and sample comparisons

## Project Structure

```
snRNA-seq-Pipeline/
├── launch.sh                    # Main launcher script
├── README.md                    # This file
├── LICENSE                      # License information
├── requirements.txt             # R package dependencies
├── .gitignore                   # Git ignore rules
├── docker-compose.yml           # Docker services configuration
├── Dockerfile                   # Container build configuration
├── scripts/                     # All executable scripts
│   ├── run_shiny_app.R         # Shiny app launcher
│   ├── run_pipeline_terminal.R  # Command-line interface
│   ├── seurat_pipeline.R        # Original pipeline script
│   ├── deployment/              # Deployment scripts
│   │   ├── deploy.sh           # Docker deployment
│   │   ├── launch_app.sh       # Shiny app launcher
│   │   ├── quick_start.sh      # Setup script
│   │   └── nginx.conf          # NGINX configuration
│   └── utils/                   # Utility scripts
│       └── test_shiny.R        # Testing utilities
├── src/                         # Source code
│   ├── core/                    # Core pipeline modules
│   ├── utils/                   # Utility functions
│   └── visualization/           # Plotting functions
├── shiny_app/                   # Shiny web application
│   ├── app.R                   # Main Shiny app
│   ├── css/                    # Styling
│   ├── js/                     # JavaScript
│   ├── www/                    # Web assets
│   └── docs/                   # App documentation
├── docs/                        # Documentation
│   ├── user_guides/            # User guides
│   │   ├── SHINY_APP_GUIDE.md  # Shiny app guide
│   │   └── DEPLOYMENT_GUIDE.md # Deployment guide
│   ├── vignettes/              # Tutorials
│   └── api/                    # API documentation
├── config/                      # Configuration files
├── data/                        # Data directory
├── results/                     # Output directory
├── logs/                        # Log files
└── temp/                        # Temporary files
```

## Quick Start

### 🚀 One-Command Setup (Recommended)

```bash
# Clone and setup in one go
git clone <repository-url>
cd snRNA-seq-Pipeline
./launch.sh setup

# Launch the Shiny app (opens in browser automatically)
./launch.sh shiny
```

### 🎯 Super Quick Start

If you just want to see the app immediately:

```bash
# Setup and launch in one go
git clone <repository-url>
cd snRNA-seq-Pipeline
./launch.sh setup && ./launch.sh shiny
```

### Option 1: Docker Deployment

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd snRNA-seq-Pipeline
   ```

2. **Deploy using the automated script**:
   ```bash
   ./launch.sh deploy
   ```

3. **Access the application**:
   - Web interface: http://localhost:3838
   - Terminal interface: Use `./launch.sh terminal` to see available commands

### Option 2: Local Installation

1. **Run the quick setup**:
   ```bash
   ./launch.sh setup
   ```

2. **Start the Shiny app**:
   ```bash
   ./launch.sh shiny
   ```

3. **Or run terminal pipeline**:
   ```bash
   ./launch.sh terminal
   ```

## Deployment Options

### Development Mode
```bash
./launch.sh deploy
```

### Production Mode
```bash
./launch.sh deploy
```

### Custom Port
```bash
# Use environment variable for custom port
SHINY_PORT=8080 ./launch.sh deploy
```

## Usage

### Web Interface

1. **Data Input**: Upload your H5 or RDS files
2. **Configure Parameters**: Set QC, processing, and clustering parameters
3. **Run Analysis**: Execute the pipeline with one click
4. **View Results**: Explore plots, tables, and downloadable files

### Command Line Interface

#### Basic Usage
```bash
# Run with H5 file
Rscript scripts/run_pipeline_terminal.R --h5_input data/sample.h5 --project_name MyProject

# Run with RDS file
Rscript scripts/run_pipeline_terminal.R --rds_input data/sample.rds --project_name MyProject

# Run with custom parameters
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name MyProject \
  --min_features 300 \
  --n_variable_features 3000 \
  --clustering_resolution 0.8 \
  --find_markers \
  --verbose
```

#### Available Parameters

**Input Options**:
- `--rds_input`: Path to Seurat RDS file
- `--h5_input`: Path to H5 data file (10X Genomics)
- `--project_name`: Project name (default: SeuratProject)
- `--working_dir`: Working directory (default: current directory)
- `--soupor_cell_doublet_input`: SouporCell output for doublet detection

**Quality Control**:
- `--min_features`: Minimum features per cell (default: 200)
- `--min_counts`: Minimum counts per cell (default: 1000)
- `--max_features`: Maximum features per cell (default: 6000)
- `--max_counts`: Maximum counts per cell (default: 25000)
- `--max_mt_percent`: Maximum mitochondrial percentage (default: 20.0)
- `--min_cells`: Minimum cells per gene (default: 3)

**Processing**:
- `--n_variable_features`: Number of variable features (default: 2000)
- `--normalization_method`: Normalization method (LogNormalize|RC|CLR, default: CLR)
- `--scaling_method`: Scaling method (negbinom|linear, default: negbinom)
- `--pca_dimensions`: PCA dimensions (default: 15)
- `--scale_factor`: Scale factor (default: 10000)

**Clustering**:
- `--clustering_resolution`: Clustering resolution (default: 0.5)
- `--clustering_algorithm`: Algorithm (louvain|multilevel|leiden|slm, default: leiden)
- `--min_cluster_size`: Minimum cluster size (default: 10)
- `--umap_n_neighbors`: UMAP neighbors (default: 30)
- `--umap_min_dist`: UMAP minimum distance (default: 0.3)

### Multi-Sample Analysis

The pipeline now supports processing multiple samples simultaneously with parallel processing and sample-specific configurations.

#### Web Interface - Multi-Sample

1. **Data Input**: Upload multiple H5 or RDS files
2. **Parallel Processing**: Enable parallel processing and set number of cores
3. **Sample Configurations**: Optionally upload sample-specific configuration file
4. **Run Analysis**: Execute the pipeline for all samples
5. **View Results**: Explore individual sample results and comparisons

#### Command Line - Multi-Sample

```bash
# Run multiple RDS files
Rscript scripts/run_multi_sample_pipeline.R \
  --rds_inputs data/sample1.rds data/sample2.rds data/sample3.rds \
  --project_name MultiSampleProject \
  --parallel \
  --n_cores 4

# Run multiple H5 files with sample configurations
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample1.h5 data/sample2.h5 \
  --project_name MultiSampleProject \
  --sample_configs config/sample_configs.yaml \
  --parallel \
  --n_cores 2

# Run with SouporCell doublet detection
Rscript scripts/run_multi_sample_pipeline.R \
  --rds_inputs data/sample1.rds data/sample2.rds \
  --soupor_cell_doublet_inputs data/soupor1.tsv data/soupor2.tsv \
  --project_name MultiSampleProject \
  --parallel
```

#### Sample-Specific Configuration

Create a YAML file to customize parameters for individual samples:

```yaml
samples:
  sample1:
    qc:
      min_features: 300
      max_mt_percent: 15
    clustering:
      resolution: 0.8
      algorithm: "louvain"
  
  sample2:
    qc:
      min_counts: 1500
    processing:
      n_variable_features: 2500
    clustering:
      resolution: 0.3
```

#### Multi-Sample Parameters

**Input Options**:
- `--rds_inputs`: Paths to multiple Seurat RDS files
- `--h5_inputs`: Paths to multiple H5 data files
- `--soupor_cell_doublet_inputs`: Paths to multiple SouporCell output files
- `--sample_configs`: Path to sample-specific configuration file

**Parallel Processing**:
- `--parallel`: Enable parallel processing (default: TRUE)
- `--n_cores`: Number of cores to use (default: auto-detect)

#### Output Structure

Multi-sample analysis creates an organized output structure:

```
ProjectName_outputs/
├── individual_samples/
│   ├── sample1/
│   │   ├── sample1_processed.rds
│   │   ├── sample1_UMAP.png
│   │   └── ...
│   ├── sample2/
│   │   ├── sample2_processed.rds
│   │   ├── sample2_UMAP.png
│   │   └── ...
│   └── ...
├── comparisons/
│   ├── ProjectName_cell_counts_comparison.png
│   ├── ProjectName_gene_counts_comparison.png
│   └── ProjectName_cluster_counts_comparison.png
├── combined_analysis/
├── reports/
└── ProjectName_multi_sample_summary.csv
```

**Analysis**:
- `--find_markers`: Find cluster markers
- `--save_intermediate`: Save intermediate results
- `--random_seed`: Random seed (default: 42)
- `--verbose`: Enable verbose output
- `--output_format`: Plot format (png|pdf|svg, default: png)

## Server Deployment

### Docker Deployment

The pipeline is containerized for easy deployment on any server:

```bash
# Build and start
docker-compose up -d

# View logs
docker-compose logs -f

# Stop services
docker-compose down
```

### Production Deployment

For production servers, use the production profile:

```bash
# Deploy with NGINX reverse proxy
docker-compose --profile production up -d

# Configure SSL certificates in nginx.conf
# Update domain name in nginx.conf
```

### Manual Server Setup

1. **Install dependencies**:
   ```bash
   sudo apt update
   sudo apt install r-base r-base-dev
   sudo apt install gdebi-core
   wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb
   sudo gdebi shiny-server-1.5.20.1002-amd64.deb
   ```

2. **Configure Shiny Server**:
   ```bash
   sudo cp shiny_app /srv/shiny-server/snrna-pipeline
   sudo systemctl restart shiny-server
   ```

3. **Set up reverse proxy** (optional):
   ```bash
   sudo apt install nginx
   sudo cp nginx.conf /etc/nginx/sites-available/snrna-pipeline
   sudo ln -s /etc/nginx/sites-available/snrna-pipeline /etc/nginx/sites-enabled/
   sudo systemctl restart nginx
   ```

## Configuration

### Environment Variables

- `SHINY_HOST`: Host address (default: 0.0.0.0)
- `SHINY_PORT`: Port number (default: 3838)
- `R_MAX_MEM_SIZE`: Maximum R memory (default: 8G)

### Configuration Files

- `config/settings.yaml`: Pipeline configuration
- `nginx.conf`: NGINX reverse proxy configuration
- `docker-compose.yml`: Docker services configuration

## Output Structure

```
results/
├── MyProject/
│   ├── qc/
│   │   ├── qc_plots.pdf
│   │   └── qc_summary.txt
│   ├── processing/
│   │   ├── variable_features.pdf
│   │   └── pca_plots.pdf
│   ├── clustering/
│   │   ├── umap_clusters.pdf
│   │   └── cluster_markers.csv
│   ├── visualizations/
│   │   ├── feature_plots.pdf
│   │   └── heatmaps.pdf
│   └── MyProject_processed.rds
```

## Troubleshooting

### Common Issues

1. **Port already in use**:
   ```bash
   SHINY_PORT=8080 ./launch.sh deploy
   ```

2. **Memory issues**:
   ```bash
   export R_MAX_MEM_SIZE=16G
   ./launch.sh deploy
   ```

3. **Permission errors**:
   ```bash
   sudo chown -R $USER:$USER .
   chmod +x launch.sh
   ```

4. **Docker issues**:
   ```bash
   ./scripts/deployment/deploy.sh clean
   ./launch.sh deploy
   ```

### Logs and Debugging

```bash
# View application logs
./scripts/deployment/deploy.sh logs

# Check service status
./scripts/deployment/deploy.sh status

# View Docker logs
docker-compose logs -f snrna-pipeline
```

## Development

### Local Development

1. **Clone and setup**:
   ```bash
   git clone <repository-url>
   cd snRNA-seq-Pipeline
   Rscript setup/install_dependencies.R
   ```

2. **Run in development mode**:
   ```bash
   ./launch.sh shiny
   ```

3. **Test terminal interface**:
   ```bash
   ./launch.sh terminal
   ```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

For issues and questions:
1. Check the troubleshooting section
2. Review the logs using `./scripts/deployment/deploy.sh logs`
3. Open an issue on the repository
4. Contact the development team

## Citation

If you use this pipeline in your research, please cite:

```
snRNA-seq Pipeline v1.0
Single-nucleus RNA sequencing analysis pipeline
https://github.com/your-repo/snRNA-seq-Pipeline
```

# snRNA-seq Pipeline

A comprehensive single-nucleus RNA sequencing analysis pipeline designed for computational biologists and researchers. Built with Seurat, this tool provides both interactive web-based analysis and high-throughput command-line processing for multi-sample studies.

## Quick Start

```bash
# Clone and setup
git clone https://github.com/VibhavSetlur/scRNA-seq_Data_Validation
cd snRNA-seq-Pipeline
./launch.sh setup

# Launch web interface
./launch.sh shiny
```

**Access**: http://localhost:3838

## Overview

This pipeline streamlines single-nucleus RNA-seq analysis with:

- **Multi-sample processing** with parallel computing
- **Sample-specific configurations** for customized analysis
- **Interactive web interface** for exploratory analysis
- **Command-line tools** for batch processing
- **Comprehensive QC** and visualization
- **Publication-ready outputs**

## Analysis Workflow

1. **Data Input**: H5 files (10X Genomics) or RDS files (Seurat objects)
2. **Quality Control**: Filtering based on features, counts, and mitochondrial content
3. **Processing**: Normalization, variable feature selection, scaling
4. **Dimensionality Reduction**: PCA and UMAP
5. **Clustering**: Multiple algorithms (Leiden, Louvain, etc.)
6. **Marker Identification**: Differential expression analysis
7. **Visualization**: Automated plot generation

## Input Data

### Supported Formats
- **H5 files**: 10X Genomics Cell Ranger output
- **RDS files**: Pre-processed Seurat objects
- **Multiple samples**: Parallel processing with sample-specific parameters

### Data Requirements
- Single-nucleus RNA-seq count matrices
- Gene expression data with cell barcodes
- Optional: SouporCell output for doublet detection

## Key Features

### Multi-Sample Analysis
Process multiple samples simultaneously with:
- Parallel computing for efficiency
- Sample-specific parameter customization
- Automated sample comparisons
- Integrated results visualization

### Quality Control
Comprehensive filtering options:
- Feature count thresholds
- UMI count filtering
- Mitochondrial content filtering
- Cell and gene quality metrics

### Advanced Analysis
- Multiple clustering algorithms
- Differential expression analysis
- Dimensionality reduction
- Marker gene identification

## Usage Options

### Web Interface (Recommended for Exploration)
```bash
./launch.sh shiny
```
- Interactive parameter adjustment
- Real-time visualization
- Sample preview and validation
- Progress monitoring

### Command Line (Recommended for Batch Processing)
```bash
# Single sample
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name MyProject

# Multiple samples
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample1.h5 data/sample2.h5 \
  --project_name MultiProject \
  --parallel \
  --n_cores 2
```

### Docker Deployment
```bash
./launch.sh deploy
```

## Output Structure

```
ProjectName_outputs/
├── individual_samples/
│   └── sample_name/
│       ├── QC plots and summaries
│       ├── Processing results
│       ├── Clustering analysis
│       └── Marker gene tables
├── comparisons/
│   ├── Sample comparisons
│   └── Cross-sample analyses
├── combined_analysis/
└── reports/
```

## Configuration

### Sample-Specific Parameters
Create YAML files to customize analysis for individual samples:

```yaml
samples:
  Control_1:
    qc:
      min_features: 200
      max_mt_percent: 20
    clustering:
      resolution: 0.5
      algorithm: "leiden"
```

See `docs/sample_configuration.md` for detailed examples.

## Installation

### System Requirements
- R 4.0+ with Bioconductor
- 8GB+ RAM recommended
- Multi-core CPU for parallel processing

### Quick Setup
```bash
./launch.sh setup
```

### Manual Installation
See `docs/installation.md` for detailed instructions.

## Documentation

- **[Installation Guide](docs/installation.md)**: Detailed setup instructions
- **[User Guide](docs/user_guide.md)**: Step-by-step analysis workflow
- **[Sample Configuration](docs/sample_configuration.md)**: Customizing parameters
- **[Advanced Usage](docs/advanced_usage.md)**: Advanced features and options
- **[Troubleshooting](docs/troubleshooting.md)**: Common issues and solutions
- **[API Reference](docs/api_reference.md)**: Function documentation

## Troubleshooting

### Common Issues
1. **Memory errors**: Increase `R_MAX_MEM_SIZE`
2. **Port conflicts**: Use `SHINY_PORT=8080 ./launch.sh shiny`
3. **Parallel processing**: Reduce `--n_cores` parameter

### Getting Help
- Check `docs/troubleshooting.md`
- Review logs: `./scripts/deployment/deploy.sh logs`
- Open an issue on GitHub

## Performance

### Recommended Settings
- **Small datasets** (<10K cells): 2-4 cores
- **Medium datasets** (10K-100K cells): 4-8 cores  
- **Large datasets** (>100K cells): 8+ cores, 16GB+ RAM

### Memory Usage
- ~2GB per 10K cells
- Scale linearly with dataset size
- Adjust `R_MAX_MEM_SIZE` as needed

## Contributing

We welcome contributions from the research community:
1. Fork the repository
2. Create a feature branch
3. Test thoroughly with your data
4. Submit a pull request

## Citation

If you use this pipeline in your research, please cite:

```
snRNA-seq Pipeline v1.0
Single-nucleus RNA sequencing analysis pipeline
https://github.com/VibhavSetlur/scRNA-seq_Data_Validation
```

## Support

- **Documentation**: Check the `docs/` folder
- **Issues**: GitHub issue tracker
- **Questions**: Open a discussion on GitHub

## License

MIT License - see [LICENSE](LICENSE) for details.

---

**For detailed documentation, examples, and advanced usage, see the `docs/` folder.**

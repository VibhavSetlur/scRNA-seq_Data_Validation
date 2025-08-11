# Troubleshooting Guide

This guide provides solutions to common issues encountered when using the snRNA-seq pipeline.

## Installation Issues

### R Package Installation Fails

**Problem**: R packages fail to install during setup.

**Solutions**:

1. **Update R and packages**:
   ```bash
   # Update R packages
   Rscript -e "update.packages(ask = FALSE)"
   
   # Reinstall dependencies
   Rscript setup/install_dependencies.R
   ```

2. **Install system dependencies** (Ubuntu/Debian):
   ```bash
   sudo apt update
   sudo apt install r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev
   ```

3. **Use conda environment**:
   ```bash
   conda create -n snrna r-base r-essentials
   conda activate snrna
   Rscript setup/install_dependencies.R
   ```

### Permission Errors

**Problem**: Permission denied when running scripts.

**Solutions**:

1. **Fix file permissions**:
   ```bash
   chmod +x launch.sh
   chmod +x scripts/*.R
   chmod +x scripts/deployment/*.sh
   ```

2. **Check ownership**:
   ```bash
   sudo chown -R $USER:$USER .
   ```

### Docker Installation Issues

**Problem**: Docker container fails to build or run.

**Solutions**:

1. **Clean Docker cache**:
   ```bash
   docker system prune -a
   docker volume prune
   ```

2. **Rebuild container**:
   ```bash
   ./scripts/deployment/deploy.sh clean
   ./launch.sh deploy
   ```

3. **Check Docker daemon**:
   ```bash
   sudo systemctl status docker
   sudo systemctl start docker
   ```

## Runtime Issues

### Memory Errors

**Problem**: "Cannot allocate vector of size X" or similar memory errors.

**Solutions**:

1. **Increase R memory limit**:
   ```bash
   export R_MAX_MEM_SIZE=16G
   ./launch.sh shiny
   ```

2. **Use swap space**:
   ```bash
   # Create swap file
   sudo fallocate -l 32G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile
   ```

3. **Reduce data size**:
   ```bash
   # Use more stringent QC filters
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --min_features 500 \
     --min_counts 2000 \
     --max_features 4000
   ```

### Port Conflicts

**Problem**: "Address already in use" when starting Shiny app.

**Solutions**:

1. **Use different port**:
   ```bash
   SHINY_PORT=8080 ./launch.sh shiny
   ```

2. **Kill existing processes**:
   ```bash
   # Find processes using port 3838
   lsof -ti:3838
   
   # Kill process
   kill -9 $(lsof -ti:3838)
   ```

3. **Check for existing Shiny processes**:
   ```bash
   ps aux | grep shiny
   pkill -f shiny
   ```

### File Not Found Errors

**Problem**: "File not found" when specifying input files.

**Solutions**:

1. **Check file paths**:
   ```bash
   # Verify file exists
   ls -la data/sample.h5
   
   # Use absolute paths
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input /full/path/to/sample.h5
   ```

2. **Check file format**:
   ```bash
   # Verify H5 file
   h5ls data/sample.h5
   
   # Check RDS file
   Rscript -e "readRDS('data/sample.rds')"
   ```

## Data Quality Issues

### Poor Clustering Results

**Problem**: Clusters are not well-separated or meaningful.

**Solutions**:

1. **Adjust QC parameters**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --min_features 300 \
     --max_mt_percent 15 \
     --min_counts 1500
   ```

2. **Increase variable features**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --n_variable_features 3000 \
     --pca_dimensions 25
   ```

3. **Try different clustering parameters**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --clustering_resolution 0.8 \
     --clustering_algorithm louvain
   ```

### High Mitochondrial Content

**Problem**: Most cells have high mitochondrial percentage.

**Solutions**:

1. **Check data quality**:
   ```r
   # Load data and check mitochondrial genes
   library(Seurat)
   seurat_obj <- readRDS("data/sample.rds")
   
   # Check mitochondrial gene names
   mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
   print(mt_genes)
   ```

2. **Adjust mitochondrial threshold**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --max_mt_percent 30
   ```

3. **Check species-specific mitochondrial genes**:
   ```r
   # For mouse data
   mt_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE)
   ```

### Low Cell Count After QC

**Problem**: Too many cells filtered out during quality control.

**Solutions**:

1. **Relax QC thresholds**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --min_features 100 \
     --min_counts 500 \
     --max_mt_percent 25
   ```

2. **Check data distribution**:
   ```r
   # Analyze QC metrics
   library(Seurat)
   seurat_obj <- readRDS("data/sample.rds")
   
   # Plot QC distributions
   VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
   ```

## Multi-Sample Analysis Issues

### Parallel Processing Errors

**Problem**: Multi-sample analysis fails with parallel processing.

**Solutions**:

1. **Disable parallel processing**:
   ```bash
   Rscript scripts/run_multi_sample_pipeline.R \
     --h5_inputs data/sample1.h5 data/sample2.h5 \
     --parallel FALSE
   ```

2. **Reduce number of cores**:
   ```bash
   Rscript scripts/run_multi_sample_pipeline.R \
     --h5_inputs data/sample1.h5 data/sample2.h5 \
     --n_cores 2
   ```

3. **Check system resources**:
   ```bash
   # Check available cores
   nproc
   
   # Check memory
   free -h
   ```

### Sample Configuration Errors

**Problem**: YAML configuration file causes errors.

**Solutions**:

1. **Validate YAML syntax**:
   ```bash
   # Check YAML syntax
   python3 -c "import yaml; yaml.safe_load(open('config/sample_configs.yaml'))"
   ```

2. **Use minimal configuration**:
   ```yaml
   samples:
     sample1:
       qc:
         min_features: 200
       clustering:
         resolution: 0.5
   ```

3. **Check parameter names**:
   - Ensure parameter names match exactly
   - Use lowercase for parameter names
   - Check for typos in section names

## Shiny App Issues

### App Won't Start

**Problem**: Shiny app fails to launch.

**Solutions**:

1. **Check R packages**:
   ```r
   # Verify required packages
   required_packages <- c("shiny", "shinydashboard", "DT", "plotly")
   for(pkg in required_packages) {
     if(!require(pkg, character.only = TRUE)) {
       install.packages(pkg)
     }
   }
   ```

2. **Check port availability**:
   ```bash
   # Test port
   netstat -tuln | grep 3838
   
   # Use different port
   SHINY_PORT=8080 ./launch.sh shiny
   ```

3. **Check file permissions**:
   ```bash
   ls -la shiny_app/
   chmod -R 755 shiny_app/
   ```

### App Crashes During Analysis

**Problem**: Shiny app crashes when running analysis.

**Solutions**:

1. **Increase memory**:
   ```bash
   export R_MAX_MEM_SIZE=16G
   ./launch.sh shiny
   ```

2. **Check logs**:
   ```bash
   # View Shiny logs
   tail -f logs/shiny.log
   
   # Check system logs
   journalctl -u shiny-server -f
   ```

3. **Use command line for large datasets**:
   ```bash
   # For large datasets, use command line instead
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/large_sample.h5 \
     --project_name LargeAnalysis
   ```

## Performance Issues

### Slow Processing

**Problem**: Analysis takes too long to complete.

**Solutions**:

1. **Use parallel processing**:
   ```bash
   Rscript scripts/run_multi_sample_pipeline.R \
     --h5_inputs data/*.h5 \
     --parallel \
     --n_cores 8
   ```

2. **Reduce data size**:
   ```bash
   # Use more stringent QC
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --min_features 500 \
     --n_variable_features 1500
   ```

3. **Use SSD storage**:
   ```bash
   # Move data to SSD if available
   cp data/sample.h5 /mnt/ssd/
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input /mnt/ssd/sample.h5
   ```

### High Memory Usage

**Problem**: Analysis uses excessive memory.

**Solutions**:

1. **Monitor memory usage**:
   ```bash
   # Monitor in real-time
   htop
   
   # Check R memory
   Rscript -e "gc(); print(memory.size())"
   ```

2. **Use memory-efficient settings**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --save_intermediate FALSE
   ```

3. **Process in chunks**:
   ```bash
   # Split large dataset
   # Process each chunk separately
   for chunk in data/chunk_*.h5; do
     Rscript scripts/run_pipeline_terminal.R \
       --h5_input "$chunk" \
       --project_name "chunk_$(basename $chunk .h5)"
   done
   ```

## Getting Help

### Debugging Steps

1. **Check logs**:
   ```bash
   # View application logs
   tail -f logs/pipeline.log
   
   # Check system logs
   journalctl -xe
   ```

2. **Enable verbose output**:
   ```bash
   Rscript scripts/run_pipeline_terminal.R \
     --h5_input data/sample.h5 \
     --verbose
   ```

3. **Test with sample data**:
   ```bash
   # Use provided test data
   ./scripts/utils/run_complete_test.sh
   ```

### Reporting Issues

When reporting issues, please include:

1. **System information**:
   ```bash
   # OS and R version
   uname -a
   R --version
   
   # Package versions
   Rscript -e "sessionInfo()"
   ```

2. **Error messages**: Copy the complete error message
3. **Command used**: The exact command that failed
4. **Data information**: File size, format, and sample information
5. **Expected vs actual behavior**: What you expected vs what happened

### Additional Resources

- [R Documentation](https://www.r-project.org/help.html)
- [Seurat Documentation](https://satijalab.org/seurat/)
- [Shiny Documentation](https://shiny.rstudio.com/)
- [GitHub Issues](https://github.com/your-repo/snRNA-seq-Pipeline/issues)

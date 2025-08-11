# Troubleshooting Guide

This guide helps you resolve common issues when using the snRNA-seq Pipeline.

## Common Issues and Solutions

### Application Won't Start

**Problem**: The Shiny app fails to start or shows an error.

**Solutions**:
1. **Check R packages**: Ensure all required packages are installed
   ```r
   Rscript setup/install_dependencies.R
   ```

2. **Check port availability**: The default port 3838 might be in use
   ```bash
   # Use a different port
   SHINY_PORT=8080 Rscript run_shiny_app.R
   ```

3. **Check file permissions**: Ensure scripts are executable
   ```bash
   chmod +x run_shiny_app.R
   chmod +x run_pipeline_terminal.R
   ```

### File Upload Issues

**Problem**: Files won't upload or show errors.

**Solutions**:
1. **Check file format**: Ensure files are in the correct format
   - H5 files: Must be 10X Genomics format
   - RDS files: Must be Seurat objects

2. **Check file size**: Large files may take time to upload
   - Wait for upload to complete
   - Check browser console for errors

3. **Check file path**: Ensure file paths don't contain special characters

### Pipeline Execution Errors

**Problem**: Pipeline fails to run or produces errors.

**Solutions**:
1. **Check input data**: Ensure data is properly formatted
   - H5 files should contain gene expression matrices
   - RDS files should be valid Seurat objects

2. **Adjust parameters**: Try more conservative settings
   - Lower minimum features/counts
   - Reduce number of variable features
   - Use fewer PCA dimensions

3. **Check memory**: Large datasets may require more memory
   ```bash
   export R_MAX_MEM_SIZE=16G
   ```

### Quality Control Issues

**Problem**: Too many cells are filtered out during QC.

**Solutions**:
1. **Relax QC thresholds**:
   - Lower minimum features/counts
   - Increase maximum features/counts
   - Increase maximum mitochondrial percentage

2. **Check data quality**: Examine QC plots to understand filtering
   - Look for unusual distributions
   - Consider if thresholds are appropriate for your data

3. **Data-specific adjustments**:
   - Different tissues may have different QC requirements
   - Adjust thresholds based on your specific dataset

### Clustering Problems

**Problem**: Clustering produces too many or too few clusters.

**Solutions**:
1. **Adjust clustering resolution**:
   - Lower resolution = fewer clusters
   - Higher resolution = more clusters

2. **Try different algorithms**:
   - leiden: Good for most datasets
   - louvain: Alternative option
   - multilevel: For specific cases

3. **Check data preprocessing**:
   - Ensure proper normalization
   - Verify variable feature selection
   - Check PCA dimensions

### Memory Issues

**Problem**: Application runs out of memory or crashes.

**Solutions**:
1. **Increase memory limit**:
   ```bash
   export R_MAX_MEM_SIZE=16G
   ```

2. **Reduce dataset size**:
   - Subsample cells if possible
   - Use fewer variable features
   - Reduce PCA dimensions

3. **Use terminal interface**: For large datasets, use the command-line interface
   ```bash
   Rscript run_pipeline_terminal.R --help
   ```

### Performance Issues

**Problem**: Analysis is very slow or unresponsive.

**Solutions**:
1. **Optimize parameters**:
   - Reduce number of variable features
   - Use fewer PCA dimensions
   - Lower clustering resolution

2. **Use parallel processing**:
   - Ensure multiple CPU cores are available
   - Check if parallel processing is enabled

3. **Consider data size**:
   - Very large datasets may require specialized hardware
   - Consider using cloud computing resources

### Output Issues

**Problem**: Output files are missing or incomplete.

**Solutions**:
1. **Check working directory**: Ensure write permissions
   ```bash
   ls -la results/
   ```

2. **Check disk space**: Ensure sufficient storage
   ```bash
   df -h
   ```

3. **Check pipeline logs**: Look for error messages
   ```bash
   ./deploy.sh logs
   ```

## Error Messages

### "Cannot open connection"
- Check if required files exist
- Verify file paths are correct
- Ensure proper permissions

### "Object not found"
- Check if required packages are installed
- Verify function names are correct
- Ensure data objects are properly loaded

### "Memory allocation failed"
- Increase memory limit
- Reduce dataset size
- Use more conservative parameters

### "Timeout error"
- Increase timeout settings
- Use smaller datasets
- Check system resources

## Getting Help

### Before Contacting Support

1. **Check logs**: Review error messages and logs
2. **Test with sample data**: Try with provided test data
3. **Document the issue**: Note exact steps and error messages
4. **Check system requirements**: Ensure your system meets requirements

### Contact Information

- **GitHub Issues**: Report bugs and request features
- **Documentation**: Check this guide and other documentation
- **Community**: Ask questions in discussions or forums

### Providing Information

When reporting issues, include:
- R version and platform
- Package versions
- Error messages
- Steps to reproduce
- Sample data (if possible)
- System specifications

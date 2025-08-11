#!/bin/bash

# Comprehensive Multi-Sample Pipeline Test Script
# This script tests the complete multi-sample functionality

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Test directory
TEST_DIR="test_multi_sample_complete"
PROJECT_NAME="TestMultiSampleComplete"

# Cleanup function
cleanup() {
    if [ -d "$TEST_DIR" ]; then
        print_status "Cleaning up test directory..."
        rm -rf "$TEST_DIR"
    fi
}

# Setup test environment
setup_test_environment() {
    print_status "Setting up test environment..."
    
    # Create test directory
    mkdir -p "$TEST_DIR"
    mkdir -p "$TEST_DIR/data"
    mkdir -p "$TEST_DIR/config"
    mkdir -p "$TEST_DIR/results"
    
    print_success "Test environment created"
}

# Generate test data
generate_test_data() {
    print_status "Generating test data..."
    
    # Create R script to generate test Seurat objects
    cat > "$TEST_DIR/generate_test_data.R" << 'EOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# Function to create mock Seurat object
create_mock_seurat <- function(sample_name, n_cells, n_genes) {
  set.seed(42)
  
  # Create random count matrix
  counts <- matrix(
    rpois(n_cells * n_genes, lambda = 2),
    nrow = n_genes,
    ncol = n_cells
  )
  
  # Create gene and cell names
  gene_names <- paste0("Gene_", 1:n_genes)
  cell_names <- paste0("Cell_", 1:n_cells)
  
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name
  )
  
  # Add metadata
  seurat_obj$sample_id <- sample_name
  seurat_obj$condition <- ifelse(grepl("Control", sample_name), "Control", "Treatment")
  
  return(seurat_obj)
}

# Create test samples
samples <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
n_cells_per_sample <- c(150, 200, 180, 220)
n_genes_per_sample <- c(2500, 3000, 2800, 3200)

for (i in seq_along(samples)) {
  sample_name <- samples[i]
  n_cells <- n_cells_per_sample[i]
  n_genes <- n_genes_per_sample[i]
  
  print(paste("Creating", sample_name, "with", n_cells, "cells and", n_genes, "genes"))
  
  seurat_obj <- create_mock_seurat(sample_name, n_cells, n_genes)
  
  # Save to file
  output_file <- file.path("data", paste0(sample_name, ".rds"))
  saveRDS(seurat_obj, file = output_file)
  
  print(paste("Saved", sample_name, "to", output_file))
}

print("Test data generation completed!")
EOF
    
    # Run the R script
    cd "$TEST_DIR"
    Rscript generate_test_data.R
    cd ..
    
    print_success "Test data generated"
}

# Create sample configuration
create_sample_config() {
    print_status "Creating sample configuration..."
    
    cat > "$TEST_DIR/config/sample_configs.yaml" << 'EOF'
# Sample-specific configuration for testing
samples:
  Control_1:
    qc:
      min_features: 200
      min_counts: 800
      max_mt_percent: 18
    clustering:
      resolution: 0.4
      algorithm: "leiden"
  
  Control_2:
    qc:
      min_features: 250
      min_counts: 1000
      max_mt_percent: 20
    clustering:
      resolution: 0.5
      algorithm: "louvain"
  
  Treatment_1:
    qc:
      min_features: 180
      min_counts: 900
      max_mt_percent: 22
    processing:
      n_variable_features: 2500
    clustering:
      resolution: 0.6
      algorithm: "leiden"
  
  Treatment_2:
    qc:
      min_features: 220
      min_counts: 1100
      max_mt_percent: 25
    processing:
      n_variable_features: 2800
    clustering:
      resolution: 0.7
      algorithm: "louvain"
EOF
    
    print_success "Sample configuration created"
}

# Test 1: Basic multi-sample pipeline
test_basic_multi_sample() {
    print_status "Test 1: Basic multi-sample pipeline..."
    
    cd "$TEST_DIR"
    
    # Run multi-sample pipeline
    Rscript ../scripts/run_multi_sample_pipeline.R \
        --rds_inputs data/Control_1.rds data/Control_2.rds data/Treatment_1.rds data/Treatment_2.rds \
        --project_name "$PROJECT_NAME" \
        --working_dir results \
        --config_file ../config/settings.yaml \
        --parallel TRUE \
        --n_cores 2 \
        --find_markers FALSE \
        --verbose TRUE
    
    # Check results
    if [ -d "results/${PROJECT_NAME}_outputs" ]; then
        print_success "Basic multi-sample pipeline completed"
    else
        print_error "Basic multi-sample pipeline failed"
        exit 1
    fi
    
    cd ..
}

# Test 2: Multi-sample with sample configurations
test_sample_configs() {
    print_status "Test 2: Multi-sample with sample configurations..."
    
    cd "$TEST_DIR"
    
    # Run multi-sample pipeline with sample configs
    Rscript ../scripts/run_multi_sample_pipeline.R \
        --rds_inputs data/Control_1.rds data/Control_2.rds data/Treatment_1.rds data/Treatment_2.rds \
        --project_name "${PROJECT_NAME}_Config" \
        --working_dir results \
        --config_file ../config/settings.yaml \
        --sample_configs config/sample_configs.yaml \
        --parallel TRUE \
        --n_cores 2 \
        --find_markers FALSE \
        --verbose TRUE
    
    # Check results
    if [ -d "results/${PROJECT_NAME}_Config_outputs" ]; then
        print_success "Multi-sample with sample configurations completed"
    else
        print_error "Multi-sample with sample configurations failed"
        exit 1
    fi
    
    cd ..
}

# Test 3: Sequential processing
test_sequential_processing() {
    print_status "Test 3: Sequential processing..."
    
    cd "$TEST_DIR"
    
    # Run multi-sample pipeline sequentially
    Rscript ../scripts/run_multi_sample_pipeline.R \
        --rds_inputs data/Control_1.rds data/Control_2.rds \
        --project_name "${PROJECT_NAME}_Sequential" \
        --working_dir results \
        --config_file ../config/settings.yaml \
        --parallel FALSE \
        --find_markers FALSE \
        --verbose TRUE
    
    # Check results
    if [ -d "results/${PROJECT_NAME}_Sequential_outputs" ]; then
        print_success "Sequential processing completed"
    else
        print_error "Sequential processing failed"
        exit 1
    fi
    
    cd ..
}

# Test 4: Error handling
test_error_handling() {
    print_status "Test 4: Error handling..."
    
    cd "$TEST_DIR"
    
    # Test with non-existent files
    if Rscript ../scripts/run_multi_sample_pipeline.R \
        --rds_inputs nonexistent1.rds nonexistent2.rds \
        --project_name "${PROJECT_NAME}_Error" \
        --working_dir results \
        --config_file ../config/settings.yaml \
        --parallel FALSE \
        --verbose TRUE 2>&1 | grep -q "Error"; then
        print_success "Error handling works correctly"
    else
        print_error "Error handling test failed"
        exit 1
    fi
    
    cd ..
}

# Test 5: Output structure validation
test_output_structure() {
    print_status "Test 5: Output structure validation..."
    
    cd "$TEST_DIR"
    
    # Check output structure
    output_dir="results/${PROJECT_NAME}_outputs"
    
    if [ ! -d "$output_dir" ]; then
        print_error "Output directory not found"
        exit 1
    fi
    
    # Check required subdirectories
    required_dirs=("individual_samples" "comparisons" "combined_analysis" "reports")
    for dir in "${required_dirs[@]}"; do
        if [ ! -d "$output_dir/$dir" ]; then
            print_error "Required subdirectory $dir not found"
            exit 1
        fi
    done
    
    # Check individual sample directories
    sample_dirs=("Control_1" "Control_2" "Treatment_1" "Treatment_2")
    for sample in "${sample_dirs[@]}"; do
        if [ ! -d "$output_dir/individual_samples/$sample" ]; then
            print_error "Sample directory $sample not found"
            exit 1
        fi
    done
    
    # Check summary file
    if [ ! -f "$output_dir/${PROJECT_NAME}_multi_sample_summary.csv" ]; then
        print_error "Summary file not found"
        exit 1
    fi
    
    print_success "Output structure validation passed"
    
    cd ..
}

# Test 6: Shiny app integration
test_shiny_integration() {
    print_status "Test 6: Shiny app integration..."
    
    # Test that Shiny app can be loaded
    if Rscript -e "
    suppressPackageStartupMessages({
      library(shiny)
      library(shinydashboard)
    })
    source('shiny_app/app.R')
    cat('Shiny app loaded successfully\n')
    " 2>&1 | grep -q "Shiny app loaded successfully"; then
        print_success "Shiny app integration test passed"
    else
        print_error "Shiny app integration test failed"
        exit 1
    fi
}

# Test 7: Performance testing
test_performance() {
    print_status "Test 7: Performance testing..."
    
    cd "$TEST_DIR"
    
    # Measure execution time
    start_time=$(date +%s)
    
    Rscript ../scripts/run_multi_sample_pipeline.R \
        --rds_inputs data/Control_1.rds data/Control_2.rds \
        --project_name "${PROJECT_NAME}_Performance" \
        --working_dir results \
        --config_file ../config/settings.yaml \
        --parallel TRUE \
        --n_cores 2 \
        --find_markers FALSE \
        --verbose FALSE
    
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    
    print_status "Performance test completed in ${duration} seconds"
    
    # Check if performance is reasonable (should complete within 5 minutes)
    if [ $duration -lt 300 ]; then
        print_success "Performance test passed"
    else
        print_warning "Performance test took longer than expected (${duration}s)"
    fi
    
    cd ..
}

# Test 8: Memory usage testing
test_memory_usage() {
    print_status "Test 8: Memory usage testing..."
    
    cd "$TEST_DIR"
    
    # Run memory test
    if command -v /usr/bin/time >/dev/null 2>&1; then
        # Use /usr/bin/time to measure memory usage
        memory_output=$(/usr/bin/time -f "%M" Rscript ../scripts/run_multi_sample_pipeline.R \
            --rds_inputs data/Control_1.rds \
            --project_name "${PROJECT_NAME}_Memory" \
            --working_dir results \
            --config_file ../config/settings.yaml \
            --parallel FALSE \
            --find_markers FALSE \
            --verbose FALSE 2>&1 | tail -1)
        
        memory_mb=$((memory_output / 1024))
        print_status "Memory usage: ${memory_mb} MB"
        
        # Check if memory usage is reasonable (less than 4GB)
        if [ $memory_mb -lt 4096 ]; then
            print_success "Memory usage test passed"
        else
            print_warning "Memory usage is high (${memory_mb} MB)"
        fi
    else
        print_warning "Memory testing skipped (time command not available)"
    fi
    
    cd ..
}

# Test 9: Cross-platform compatibility
test_cross_platform() {
    print_status "Test 9: Cross-platform compatibility..."
    
    # Test R version
    r_version=$(R --version | head -1)
    print_status "R version: $r_version"
    
    # Test platform
    platform=$(R -e "cat(R.version\$platform)" 2>/dev/null)
    print_status "Platform: $platform"
    
    # Test required packages
    missing_packages=$(R -e "
    required_packages <- c('Seurat', 'future', 'future.apply', 'parallel', 'purrr', 'yaml', 'fs')
    missing <- character()
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        missing <- c(missing, pkg)
      }
    }
    if (length(missing) > 0) {
      cat(paste(missing, collapse = ', '))
    } else {
      cat('All packages available')
    }
    " 2>/dev/null)
    
    if [ "$missing_packages" = "All packages available" ]; then
        print_success "Cross-platform compatibility test passed"
    else
        print_error "Missing packages: $missing_packages"
        exit 1
    fi
}

# Test 10: Documentation and examples
test_documentation() {
    print_status "Test 10: Documentation and examples..."
    
    # Check if documentation files exist
    required_files=(
        "README.md"
        "config/sample_configs_template.yaml"
        "scripts/run_multi_sample_pipeline.R"
        "tests/test_multi_sample_pipeline.R"
    )
    
    for file in "${required_files[@]}"; do
        if [ -f "$file" ]; then
            print_success "Documentation file exists: $file"
        else
            print_error "Documentation file missing: $file"
            exit 1
        fi
    done
    
    # Check if examples are valid
    if grep -q "multi-sample" README.md; then
        print_success "Documentation contains multi-sample examples"
    else
        print_error "Documentation missing multi-sample examples"
        exit 1
    fi
}

# Main test execution
main() {
    print_status "=== Starting Comprehensive Multi-Sample Pipeline Tests ==="
    
    # Setup
    cleanup
    setup_test_environment
    generate_test_data
    create_sample_config
    
    # Run tests
    test_basic_multi_sample
    test_sample_configs
    test_sequential_processing
    test_error_handling
    test_output_structure
    test_shiny_integration
    test_performance
    test_memory_usage
    test_cross_platform
    test_documentation
    
    # Final cleanup
    cleanup
    
    print_success "=== All Multi-Sample Pipeline Tests Completed Successfully ==="
    print_status "The multi-sample functionality is working correctly!"
}

# Run main function
main "$@"

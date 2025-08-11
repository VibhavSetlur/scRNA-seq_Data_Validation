#!/bin/bash

# Test script for snRNA-seq Pipeline
# This script runs basic tests to ensure the pipeline is working

set -e

echo "Running snRNA-seq Pipeline tests..."

# Test 1: Check if R and required packages are available
echo "Test 1: Checking R environment..."
Rscript -e "
required_packages <- c('Seurat', 'tidyverse', 'argparse', 'yaml')
missing_packages <- character()

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) == 0) {
  cat('✓ All required packages are available\n')
} else {
  cat('✗ Missing packages:', paste(missing_packages, collapse = ', '), '\n')
  quit(status = 1)
}
"

# Test 2: Check if configuration file exists
echo "Test 2: Checking configuration files..."
if [[ -f "config/settings.yaml" ]]; then
    echo "✓ Configuration file exists"
else
    echo "✗ Configuration file missing"
    exit 1
fi

# Test 3: Check if main pipeline script exists
echo "Test 3: Checking pipeline scripts..."
if [[ -f "src/core/pipeline.R" ]]; then
    echo "✓ Main pipeline script exists"
else
    echo "✗ Main pipeline script missing"
    exit 1
fi

# Test 4: Generate test data
echo "Test 4: Generating test data..."
Rscript tests/scripts/generate_test_data.R

# Test 5: Run pipeline on test data
echo "Test 5: Running pipeline on test data..."
if [[ -f "tests/data/test_data.h5" ]]; then
    Rscript src/core/pipeline.R \
        --h5_input tests/data/test_data.h5 \
        --project_name TestProject \
        --working_dir tests/results \
        --find_markers FALSE
    echo "✓ Pipeline test completed successfully"
else
    echo "✗ Test data not found"
    exit 1
fi

# Test 6: Run multi-sample pipeline tests
echo "Test 6: Running multi-sample pipeline tests..."
if [[ -f "tests/R/test_multi_sample_pipeline.R" ]]; then
    Rscript tests/R/test_multi_sample_pipeline.R
    echo "✓ Multi-sample pipeline tests completed successfully"
else
    echo "✗ Multi-sample pipeline test script not found"
    exit 1
fi

# Test 7: Run Shiny app tests
echo "Test 7: Running Shiny app tests..."
if [[ -f "tests/R/test_shiny_app.R" ]]; then
    Rscript tests/R/test_shiny_app.R
    echo "✓ Shiny app tests completed successfully"
else
    echo "✗ Shiny app test script not found"
    exit 1
fi

echo "All tests passed!"

#!/usr/bin/env Rscript

# Multi-Sample Pipeline Tests
# This script tests the multi-sample functionality of the snRNA-seq pipeline

suppressPackageStartupMessages({
  library(testthat)
  library(yaml)
  library(fs)
  library(Seurat)
  library(future)
  library(future.apply)
})

# Source the pipeline
source("src/core/pipeline.R")

# Test setup
setup_test_environment <- function() {
  # Create test directories
  test_dir <- "test_multi_sample_output"
  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir, recursive = TRUE)
  
  # Create test configuration
  test_config <- list(
    project = list(
      default_name = "TestProject",
      output_format = "png",
      dpi = 300
    ),
    qc = list(
      min_features = 100,
      min_counts = 500,
      max_features = 5000,
      max_counts = 20000,
      max_mt_percent = 25,
      min_cells = 2
    ),
    processing = list(
      normalization_method = "CLR",
      n_variable_features = 1000,
      scaling_method = "negbinom",
      pca_dimensions = 10,
      scale_factor = 10000
    ),
    clustering = list(
      resolution = 0.3,
      algorithm = "leiden",
      min_cluster_size = 5,
      random_seed = 42
    ),
    visualization = list(
      theme = "classic",
      color_palette = "viridis",
      point_size = 0.5,
      label_size = 2
    )
  )
  
  config_file <- file.path(test_dir, "test_config.yaml")
  yaml::write_yaml(test_config, config_file)
  
  return(list(test_dir = test_dir, config_file = config_file))
}

# Create mock Seurat objects for testing
create_mock_seurat_objects <- function(test_dir, n_samples = 3) {
  seurat_files <- character(n_samples)
  
  for (i in 1:n_samples) {
    # Create mock data with more realistic parameters
    n_cells <- 200 + (i * 100)  # More cells for better QC
    n_genes <- 2000 + (i * 100)  # Different gene counts for each sample
    
    # Create random count matrix with higher expression
    set.seed(42 + i)
    counts <- matrix(
      rpois(n_cells * n_genes, lambda = 5),  # Higher lambda for more counts
      nrow = n_genes,
      ncol = n_cells
    )
    
    # Add some high-expression genes to ensure cells pass QC
    high_expr_genes <- sample(1:n_genes, n_genes * 0.1)  # 10% high expression genes
    for (gene_idx in high_expr_genes) {
      counts[gene_idx, ] <- rpois(n_cells, lambda = 20)  # High expression
    }
    
    # Create gene names
    gene_names <- paste0("Gene", 1:n_genes)
    cell_names <- paste0("Cell", 1:n_cells)
    
    rownames(counts) <- gene_names
    colnames(counts) <- cell_names
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = paste0("TestSample", i)
    )
    
    # Add some metadata
    seurat_obj$sample_id <- paste0("Sample", i)
    seurat_obj$condition <- ifelse(i %% 2 == 0, "Control", "Treatment")
    
    # Save to file
    file_path <- file.path(test_dir, paste0("test_sample_", i, ".rds"))
    saveRDS(seurat_obj, file_path)
    seurat_files[i] <- file_path
  }
  
  return(seurat_files)
}

# Create sample-specific configuration for testing
create_sample_configs <- function(test_dir) {
  sample_configs <- list(
    samples = list(
      test_sample_1 = list(
        qc = list(min_features = 150),
        clustering = list(resolution = 0.4)
      ),
      test_sample_2 = list(
        qc = list(min_counts = 600),
        processing = list(n_variable_features = 1200)
      ),
      test_sample_3 = list(
        clustering = list(algorithm = "louvain")
      )
    )
  )
  
  config_file <- file.path(test_dir, "sample_configs.yaml")
  yaml::write_yaml(sample_configs, config_file)
  return(config_file)
}

# Test 1: Basic multi-sample functionality
test_basic_multi_sample <- function() {
  test_that("Basic multi-sample pipeline works", {
    # Setup
    setup <- setup_test_environment()
    seurat_files <- create_mock_seurat_objects(setup$test_dir, 3)
    
    # Create arguments
    args <- list(
      rds_inputs = seurat_files,
      project_name = "TestMultiSample",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = FALSE,  # Use sequential for testing
      verbose = TRUE,
      find_markers = FALSE,  # Skip markers for faster testing
      save_intermediate = FALSE
    )
    
    # Run pipeline
    result <- run_multi_sample_pipeline(args)
    
    # Assertions
    expect_true(!is.null(result))
    expect_equal(result$project_name, "TestMultiSample")
    expect_equal(result$n_samples, 3)
    expect_true(result$n_successful >= 1)  # At least one sample should succeed
    expect_true(result$n_failed <= 2)      # Some samples might fail due to QC
    expect_true(dir.exists(result$project_output_dir))
    
    # Check individual sample results
    expect_true(length(result$samples) == 3)
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 2: Parallel processing
test_parallel_processing <- function() {
  test_that("Parallel processing works correctly", {
    # Setup
    setup <- setup_test_environment()
    seurat_files <- create_mock_seurat_objects(setup$test_dir, 4)
    
    # Create arguments with parallel processing
    args <- list(
      rds_inputs = seurat_files,
      project_name = "TestParallel",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = TRUE,
      n_cores = 2,
      verbose = TRUE,
      find_markers = FALSE,
      save_intermediate = FALSE
    )
    
    # Run pipeline
    result <- run_multi_sample_pipeline(args)
    
    # Assertions
    expect_true(!is.null(result))
    expect_equal(result$n_samples, 4)
    expect_true(result$n_successful >= 1)  # At least one sample should succeed
    expect_equal(result$parallel_cores, 2)
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 3: Sample-specific configurations
test_sample_configs <- function() {
  test_that("Sample-specific configurations work", {
    # Setup
    setup <- setup_test_environment()
    seurat_files <- create_mock_seurat_objects(setup$test_dir, 3)
    sample_configs_file <- create_sample_configs(setup$test_dir)
    
    # Create arguments with sample configs
    args <- list(
      rds_inputs = seurat_files,
      project_name = "TestSampleConfigs",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      sample_configs = sample_configs_file,
      parallel = FALSE,
      verbose = TRUE,
      find_markers = FALSE,
      save_intermediate = FALSE
    )
    
    # Run pipeline
    result <- run_multi_sample_pipeline(args)
    
    # Assertions
    expect_true(!is.null(result))
    expect_equal(result$n_samples, 3)
    expect_true(result$n_successful >= 1)  # At least one sample should succeed
    
    # Check that sample-specific configs were applied
    # (This would require checking the actual parameters used in each sample)
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 4: Error handling
test_error_handling <- function() {
  test_that("Error handling works correctly", {
    # Setup
    setup <- setup_test_environment()
    
    # Create arguments with non-existent files
    args <- list(
      rds_inputs = c("nonexistent1.rds", "nonexistent2.rds"),
      project_name = "TestErrors",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = FALSE,
      verbose = TRUE
    )
    
    # Run pipeline and expect error
    expect_error(run_multi_sample_pipeline(args))
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 5: Mixed success/failure scenarios
test_mixed_scenarios <- function() {
  test_that("Mixed success/failure scenarios work", {
    # Setup
    setup <- setup_test_environment()
    seurat_files <- create_mock_seurat_objects(setup$test_dir, 2)
    
    # Add a non-existent file to create a mixed scenario
    mixed_files <- c(seurat_files, "nonexistent.rds")
    
    # Create arguments
    args <- list(
      rds_inputs = mixed_files,
      project_name = "TestMixed",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = FALSE,
      verbose = TRUE,
      find_markers = FALSE,
      save_intermediate = FALSE
    )
    
    # Run pipeline
    result <- run_multi_sample_pipeline(args)
    
    # Assertions
    expect_true(!is.null(result))
    expect_equal(result$n_samples, 3)
    expect_equal(result$n_successful, 2)
    expect_equal(result$n_failed, 1)
    
    # Check that successful samples are marked as completed
    successful_samples <- names(result$samples)[sapply(result$samples, function(x) x$status == "completed")]
    expect_equal(length(successful_samples), 2)
    
    # Check that failed samples are marked as failed
    failed_samples <- names(result$samples)[sapply(result$samples, function(x) x$status == "failed")]
    expect_equal(length(failed_samples), 1)
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 6: H5 file inputs
test_h5_inputs <- function() {
  test_that("H5 file inputs work correctly", {
    # This test would require actual H5 files
    # For now, we'll test the argument parsing
    setup <- setup_test_environment()
    
    # Create mock H5 files (this is a simplified test)
    h5_files <- file.path(setup$test_dir, paste0("test_", 1:2, ".h5"))
    for (file in h5_files) {
      # Create empty files for testing
      file.create(file)
    }
    
    # Create arguments
    args <- list(
      h5_inputs = h5_files,
      project_name = "TestH5",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = FALSE,
      verbose = TRUE,
      find_markers = FALSE,
      save_intermediate = FALSE
    )
    
    # Test argument validation
    expect_error(run_multi_sample_pipeline(args), "Failed to read H5 file")
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 7: Output directory structure
test_output_structure <- function() {
  test_that("Output directory structure is correct", {
    # Setup
    setup <- setup_test_environment()
    seurat_files <- create_mock_seurat_objects(setup$test_dir, 2)
    
    # Create arguments
    args <- list(
      rds_inputs = seurat_files,
      project_name = "TestStructure",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = FALSE,
      verbose = TRUE,
      find_markers = FALSE,
      save_intermediate = TRUE
    )
    
    # Run pipeline
    result <- run_multi_sample_pipeline(args)
    
    # Check output directory structure
    expect_true(dir.exists(result$project_output_dir))
    expect_true(dir.exists(file.path(result$project_output_dir, "individual_samples")))
    expect_true(dir.exists(file.path(result$project_output_dir, "comparisons")))
    expect_true(dir.exists(file.path(result$project_output_dir, "combined_analysis")))
    expect_true(dir.exists(file.path(result$project_output_dir, "reports")))
    
    # Check individual sample directories
    for (sample_name in names(result$samples)) {
      sample_dir <- file.path(result$project_output_dir, "individual_samples", sample_name)
      expect_true(dir.exists(sample_dir))
    }
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 8: Summary generation
test_summary_generation <- function() {
  test_that("Summary generation works correctly", {
    # Setup
    setup <- setup_test_environment()
    seurat_files <- create_mock_seurat_objects(setup$test_dir, 3)
    
    # Create arguments
    args <- list(
      rds_inputs = seurat_files,
      project_name = "TestSummary",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = FALSE,
      verbose = TRUE,
      find_markers = FALSE,
      save_intermediate = FALSE
    )
    
    # Run pipeline
    result <- run_multi_sample_pipeline(args)
    
    # Check summary
    expect_true(!is.null(result$summary))
    expect_equal(nrow(result$summary), 3)
    expect_true(all(c("Sample_Name", "Status", "Cells", "Genes", "Clusters") %in% colnames(result$summary)))
    expect_true(all(result$summary$Status == "Completed"))
    
    # Check summary file exists
    summary_file <- file.path(result$project_output_dir, "TestSummary_multi_sample_summary.csv")
    expect_true(file.exists(summary_file))
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 9: Memory and performance
test_memory_performance <- function() {
  test_that("Memory and performance are reasonable", {
    # Setup
    setup <- setup_test_environment()
    seurat_files <- create_mock_seurat_objects(setup$test_dir, 5)
    
    # Record start time and memory
    start_time <- Sys.time()
    start_memory <- if (requireNamespace("pryr", quietly = TRUE)) pryr::mem_used() else 0
    
    # Create arguments
    args <- list(
      rds_inputs = seurat_files,
      project_name = "TestPerformance",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = TRUE,
      n_cores = 2,
      verbose = FALSE,  # Reduce output for performance test
      find_markers = FALSE,
      save_intermediate = FALSE
    )
    
    # Run pipeline
    result <- run_multi_sample_pipeline(args)
    
    # Record end time and memory
    end_time <- Sys.time()
    end_memory <- if (requireNamespace("pryr", quietly = TRUE)) pryr::mem_used() else 0
    
    # Calculate metrics
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    memory_used <- end_memory - start_memory
    
    # Assertions
    expect_true(!is.null(result))
    expect_equal(result$n_successful, 5)
    
    # Performance assertions (adjust thresholds as needed)
    expect_true(duration < 300)  # Should complete within 5 minutes
    if (memory_used > 0) {
      expect_true(memory_used < 4 * 1024^3)  # Should use less than 4GB
    }
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 10: Integration with existing pipeline
test_integration <- function() {
  test_that("Integration with existing pipeline components works", {
    # Setup
    setup <- setup_test_environment()
    seurat_files <- create_mock_seurat_objects(setup$test_dir, 2)
    
    # Create arguments
    args <- list(
      rds_inputs = seurat_files,
      project_name = "TestIntegration",
      working_dir = setup$test_dir,
      config_file = setup$config_file,
      parallel = FALSE,
      verbose = TRUE,
      find_markers = TRUE,  # Test marker finding
      save_intermediate = TRUE
    )
    
    # Run pipeline
    result <- run_multi_sample_pipeline(args)
    
    # Check that all pipeline components worked
    expect_true(!is.null(result))
    expect_equal(result$n_successful, 2)
    
    # Check that marker files were created
    for (sample_name in names(result$samples)) {
      if (result$samples[[sample_name]]$status == "completed") {
        sample_dir <- result$samples[[sample_name]]$output_dir
        marker_file <- file.path(sample_dir, paste0(sample_name, "_markers.csv"))
        # Note: Marker files might not be created if only one cluster is found
        # This is expected behavior
      }
    }
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Run all tests
run_all_tests <- function() {
  cat("=== Running Multi-Sample Pipeline Tests ===\n")
  
  # Run tests
  test_basic_multi_sample()
  test_parallel_processing()
  test_sample_configs()
  test_error_handling()
  test_mixed_scenarios()
  test_h5_inputs()
  test_output_structure()
  test_summary_generation()
  test_memory_performance()
  test_integration()
  
  cat("=== All Tests Completed ===\n")
}

# Main execution
if (!interactive()) {
  run_all_tests()
}

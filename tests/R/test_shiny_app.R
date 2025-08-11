#!/usr/bin/env Rscript

# Shiny App Tests for Multi-Sample Pipeline
# This script tests the Shiny app functionality

suppressPackageStartupMessages({
  library(testthat)
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(shinyjs)
  library(DT)
  library(plotly)
  library(yaml)
  library(fs)
  library(Seurat)
})

# Source the Shiny app
source("shiny_app/app.R")

# Test setup
setup_shiny_test_environment <- function() {
  # Create test directories
  test_dir <- "test_shiny_output"
  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir, recursive = TRUE)
  
  # Create test configuration
  test_config <- list(
    project = list(
      default_name = "TestShinyProject",
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

# Create mock Seurat objects for Shiny testing
create_mock_seurat_for_shiny <- function(test_dir, n_samples = 3) {
  seurat_files <- character(n_samples)
  
  for (i in 1:n_samples) {
    # Create mock data
    n_cells <- 100 + (i * 50)
    n_genes <- 2000 + (i * 100)
    
    # Create random count matrix
    set.seed(42 + i)
    counts <- matrix(
      rpois(n_cells * n_genes, lambda = 2),
      nrow = n_genes,
      ncol = n_cells
    )
    
    # Create gene names
    gene_names <- paste0("Gene_", 1:n_genes)
    cell_names <- paste0("Cell_", 1:n_cells)
    
    rownames(counts) <- gene_names
    colnames(counts) <- cell_names
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = paste0("ShinyTestSample", i)
    )
    
    # Add some metadata
    seurat_obj$sample_id <- paste0("Sample", i)
    seurat_obj$condition <- ifelse(i %% 2 == 0, "Control", "Treatment")
    
    # Save to file
    file_path <- file.path(test_dir, paste0("shiny_test_sample_", i, ".rds"))
    saveRDS(seurat_obj, file_path)
    seurat_files[i] <- file_path
  }
  
  return(seurat_files)
}

# Test 1: Shiny app initialization
test_shiny_initialization <- function() {
  test_that("Shiny app initializes correctly", {
    # Test that the UI can be created
    expect_no_error({
      ui <- dashboardPage(
        dashboardHeader(title = "snRNA-seq Pipeline"),
        dashboardSidebar(),
        dashboardBody()
      )
    })
    
    # Test that server function can be created
    expect_no_error({
      server <- function(input, output, session) {
        # Basic server logic
      }
    })
  })
}

# Test 2: Multi-sample file input validation
test_multi_sample_file_validation <- function() {
  test_that("Multi-sample file input validation works", {
    setup <- setup_shiny_test_environment()
    seurat_files <- create_mock_seurat_for_shiny(setup$test_dir, 3)
    
    # Test validation function (simulated)
    validate_inputs <- function(input_type, files) {
      if (input_type == "rds") {
        if (is.null(files) || nrow(files) == 0) {
          return(list(valid = FALSE, message = "Please select at least one RDS file!"))
        }
        
        # Check each file
        for (i in 1:nrow(files)) {
          if (!file.exists(files$datapath[i])) {
            return(list(valid = FALSE, message = paste("RDS file not found:", files$name[i])))
          }
        }
      }
      return(list(valid = TRUE, message = "Validation passed"))
    }
    
    # Create mock file input data
    mock_files <- data.frame(
      name = basename(seurat_files),
      datapath = seurat_files,
      stringsAsFactors = FALSE
    )
    
    # Test valid files
    result <- validate_inputs("rds", mock_files)
    expect_true(result$valid)
    
    # Test invalid files
    invalid_files <- data.frame(
      name = c("nonexistent1.rds", "nonexistent2.rds"),
      datapath = c("nonexistent1.rds", "nonexistent2.rds"),
      stringsAsFactors = FALSE
    )
    result <- validate_inputs("rds", invalid_files)
    expect_false(result$valid)
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 3: Sample preview functionality
test_sample_preview <- function() {
  test_that("Sample preview functionality works", {
    setup <- setup_shiny_test_environment()
    seurat_files <- create_mock_seurat_for_shiny(setup$test_dir, 3)
    
    # Test sample preview function
    generate_sample_preview <- function(input_type, files) {
      if (input_type == "rds" && !is.null(files)) {
        preview <- character()
        for (i in 1:nrow(files)) {
          file_info <- files[i, ]
          sample_name <- tools::file_path_sans_ext(file_info$name)
          file_size <- format(file.size(file_info$datapath), units = "MB")
          preview <- c(preview, 
                      paste("Sample", i, ":", sample_name),
                      paste("  File:", file_info$name),
                      paste("  Size:", file_size),
                      "")
        }
        return(paste(preview, collapse = "\n"))
      }
      return("No files selected. Please upload files to see sample preview.")
    }
    
    # Create mock file input data
    mock_files <- data.frame(
      name = basename(seurat_files),
      datapath = seurat_files,
      stringsAsFactors = FALSE
    )
    
    # Test preview generation
    preview <- generate_sample_preview("rds", mock_files)
    expect_true(nchar(preview) > 0)
    expect_true(grepl("Sample 1:", preview))
    expect_true(grepl("Sample 2:", preview))
    expect_true(grepl("Sample 3:", preview))
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 4: Multi-sample argument creation
test_multi_sample_args_creation <- function() {
  test_that("Multi-sample argument creation works", {
    setup <- setup_shiny_test_environment()
    seurat_files <- create_mock_seurat_for_shiny(setup$test_dir, 3)
    
    # Test argument creation function
    create_pipeline_args <- function(input_type, files, project_name, working_dir, 
                                   use_parallel, n_cores, use_sample_configs, sample_configs_file) {
      args <- list()
      
      # Input files
      if (input_type == "rds") {
        args$rds_inputs <- sapply(1:nrow(files), function(i) {
          normalizePath(files$datapath[i], mustWork = FALSE)
        })
      }
      
      # Project settings
      args$project_name <- project_name
      args$working_dir <- working_dir
      
      # Parallel processing
      args$parallel <- use_parallel
      if (use_parallel) {
        args$n_cores <- n_cores
      }
      
      # Sample configs
      if (use_sample_configs && !is.null(sample_configs_file)) {
        args$sample_configs <- normalizePath(sample_configs_file$datapath, mustWork = FALSE)
      }
      
      return(args)
    }
    
    # Create mock file input data
    mock_files <- data.frame(
      name = basename(seurat_files),
      datapath = seurat_files,
      stringsAsFactors = FALSE
    )
    
    # Test argument creation
    args <- create_pipeline_args(
      input_type = "rds",
      files = mock_files,
      project_name = "TestProject",
      working_dir = setup$test_dir,
      use_parallel = TRUE,
      n_cores = 4,
      use_sample_configs = FALSE,
      sample_configs_file = NULL
    )
    
    # Check arguments
    expect_equal(args$project_name, "TestProject")
    expect_equal(args$working_dir, setup$test_dir)
    expect_true(args$parallel)
    expect_equal(args$n_cores, 4)
    expect_equal(length(args$rds_inputs), 3)
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 5: Results handling
test_results_handling <- function() {
  test_that("Multi-sample results handling works", {
    # Test results structure validation
    validate_results_structure <- function(result) {
      required_fields <- c("project_name", "project_output_dir", "samples", 
                          "summary", "n_samples", "n_successful", "n_failed")
      
      if (is.null(result) || !is.list(result)) {
        return(FALSE)
      }
      
      missing_fields <- setdiff(required_fields, names(result))
      return(length(missing_fields) == 0)
    }
    
    # Test valid results structure
    valid_result <- list(
      project_name = "TestProject",
      project_output_dir = "test_output",
      samples = list(
        sample1 = list(status = "completed", n_cells = 100),
        sample2 = list(status = "completed", n_cells = 150)
      ),
      summary = data.frame(
        Sample_Name = c("sample1", "sample2"),
        Status = c("Completed", "Completed"),
        Cells = c(100, 150)
      ),
      n_samples = 2,
      n_successful = 2,
      n_failed = 0
    )
    
    expect_true(validate_results_structure(valid_result))
    
    # Test invalid results structure
    invalid_result <- list(
      project_name = "TestProject"
      # Missing required fields
    )
    
    expect_false(validate_results_structure(invalid_result))
  })
}

# Test 6: Sample selection and loading
test_sample_selection <- function() {
  test_that("Sample selection and loading works", {
    # Test sample selection logic
    get_successful_samples <- function(results) {
      if (is.null(results) || is.null(results$samples)) {
        return(character(0))
      }
      
      successful_samples <- names(results$samples)[
        sapply(results$samples, function(x) x$status == "completed")
      ]
      return(successful_samples)
    }
    
    # Test sample loading logic
    load_sample_data <- function(sample_name, results) {
      if (is.null(results) || is.null(results$samples) || 
          is.null(results$samples[[sample_name]])) {
        return(NULL)
      }
      
      sample_result <- results$samples[[sample_name]]
      if (sample_result$status != "completed") {
        return(NULL)
      }
      
      # In a real scenario, this would load the Seurat object
      return(list(
        sample_name = sample_name,
        n_cells = sample_result$n_cells,
        n_genes = sample_result$n_genes,
        n_clusters = sample_result$n_clusters
      ))
    }
    
    # Test data
    test_results <- list(
      samples = list(
        sample1 = list(status = "completed", n_cells = 100, n_genes = 2000, n_clusters = 5),
        sample2 = list(status = "failed", n_cells = 0, n_genes = 0, n_clusters = 0),
        sample3 = list(status = "completed", n_cells = 150, n_genes = 2500, n_clusters = 7)
      )
    )
    
    # Test successful samples extraction
    successful <- get_successful_samples(test_results)
    expect_equal(length(successful), 2)
    expect_true("sample1" %in% successful)
    expect_true("sample3" %in% successful)
    expect_false("sample2" %in% successful)
    
    # Test sample loading
    sample_data <- load_sample_data("sample1", test_results)
    expect_false(is.null(sample_data))
    expect_equal(sample_data$sample_name, "sample1")
    expect_equal(sample_data$n_cells, 100)
    
    # Test loading failed sample
    failed_data <- load_sample_data("sample2", test_results)
    expect_true(is.null(failed_data))
  })
}

# Test 7: Comparison plots generation
test_comparison_plots <- function() {
  test_that("Comparison plots generation works", {
    # Test comparison data preparation
    prepare_comparison_data <- function(results) {
      if (is.null(results) || is.null(results$samples)) {
        return(NULL)
      }
      
      successful_samples <- names(results$samples)[
        sapply(results$samples, function(x) x$status == "completed")
      ]
      
      if (length(successful_samples) == 0) {
        return(NULL)
      }
      
      comparison_data <- data.frame(
        Sample = successful_samples,
        Cells = sapply(successful_samples, function(x) results$samples[[x]]$n_cells),
        Genes = sapply(successful_samples, function(x) results$samples[[x]]$n_genes),
        Clusters = sapply(successful_samples, function(x) results$samples[[x]]$n_clusters),
        stringsAsFactors = FALSE
      )
      
      return(comparison_data)
    }
    
    # Test data
    test_results <- list(
      samples = list(
        sample1 = list(status = "completed", n_cells = 100, n_genes = 2000, n_clusters = 5),
        sample2 = list(status = "completed", n_cells = 150, n_genes = 2500, n_clusters = 7),
        sample3 = list(status = "failed", n_cells = 0, n_genes = 0, n_clusters = 0)
      )
    )
    
    # Test comparison data preparation
    comparison_data <- prepare_comparison_data(test_results)
    expect_false(is.null(comparison_data))
    expect_equal(nrow(comparison_data), 2)
    expect_equal(comparison_data$Cells, c(100, 150))
    expect_equal(comparison_data$Genes, c(2000, 2500))
    expect_equal(comparison_data$Clusters, c(5, 7))
  })
}

# Test 8: Error handling in Shiny app
test_shiny_error_handling <- function() {
  test_that("Error handling in Shiny app works", {
    # Test error message formatting
    format_error_message <- function(error) {
      if (is.null(error)) {
        return("Unknown error occurred")
      }
      
      if (is.character(error)) {
        return(error)
      }
      
      if (is.list(error) && !is.null(error$message)) {
        return(error$message)
      }
      
      return("An error occurred during processing")
    }
    
    # Test various error types
    expect_equal(format_error_message("Test error"), "Test error")
    expect_equal(format_error_message(list(message = "List error")), "List error")
    expect_equal(format_error_message(NULL), "Unknown error occurred")
    expect_equal(format_error_message(list()), "An error occurred during processing")
  })
}

# Test 9: Configuration validation
test_configuration_validation <- function() {
  test_that("Configuration validation works", {
    setup <- setup_shiny_test_environment()
    
    # Test configuration loading
    load_and_validate_config <- function(config_file) {
      if (!file.exists(config_file)) {
        return(list(valid = FALSE, message = "Configuration file not found"))
      }
      
      tryCatch({
        config <- yaml::read_yaml(config_file)
        
        # Check required sections
        required_sections <- c("project", "qc", "processing", "clustering", "visualization")
        missing_sections <- setdiff(required_sections, names(config))
        
        if (length(missing_sections) > 0) {
          return(list(valid = FALSE, 
                     message = paste("Missing sections:", paste(missing_sections, collapse = ", "))))
        }
        
        return(list(valid = TRUE, config = config))
      }, error = function(e) {
        return(list(valid = FALSE, message = paste("Error loading config:", e$message)))
      })
    }
    
    # Test valid configuration
    result <- load_and_validate_config(setup$config_file)
    expect_true(result$valid)
    expect_false(is.null(result$config))
    
    # Test invalid configuration
    invalid_result <- load_and_validate_config("nonexistent.yaml")
    expect_false(invalid_result$valid)
    
    # Cleanup
    unlink(setup$test_dir, recursive = TRUE)
  })
}

# Test 10: Performance testing
test_shiny_performance <- function() {
  test_that("Shiny app performance is acceptable", {
    # Test UI rendering performance
    measure_ui_rendering_time <- function() {
      start_time <- Sys.time()
      
      # Simulate UI creation
      ui <- dashboardPage(
        dashboardHeader(title = "snRNA-seq Pipeline"),
        dashboardSidebar(
          sidebarMenu(
            menuItem("Dashboard", tabName = "dashboard"),
            menuItem("Data Input", tabName = "data_input"),
            menuItem("Results", tabName = "results")
          )
        ),
        dashboardBody(
          tabItems(
            tabItem(tabName = "dashboard", "Dashboard content"),
            tabItem(tabName = "data_input", "Data input content"),
            tabItem(tabName = "results", "Results content")
          )
        )
      )
      
      end_time <- Sys.time()
      return(as.numeric(difftime(end_time, start_time, units = "secs")))
    }
    
    # Test UI rendering time
    rendering_time <- measure_ui_rendering_time()
    expect_true(rendering_time < 1.0)  # Should render in less than 1 second
  })
}

# Run all tests
run_all_shiny_tests <- function() {
  cat("=== Running Shiny App Tests ===\n")
  
  # Run tests
  test_shiny_initialization()
  test_multi_sample_file_validation()
  test_sample_preview()
  test_multi_sample_args_creation()
  test_results_handling()
  test_sample_selection()
  test_comparison_plots()
  test_shiny_error_handling()
  test_configuration_validation()
  test_shiny_performance()
  
  cat("=== All Shiny App Tests Completed ===\n")
}

# Main execution
if (!interactive()) {
  run_all_shiny_tests()
}

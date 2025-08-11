#!/usr/bin/env Rscript

# Multi-Sample snRNA-seq Pipeline Runner
# This script provides a command-line interface for running multi-sample analysis

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(fs)
})

# Source the pipeline
# Try multiple possible paths
pipeline_paths <- c(
  "src/core/pipeline.R",
  "../src/core/pipeline.R",
  file.path(dirname(getwd()), "src", "core", "pipeline.R")
)

pipeline_found <- FALSE
for (path in pipeline_paths) {
  if (file.exists(path)) {
    source(path)
    pipeline_found <- TRUE
    break
  }
}

if (!pipeline_found) {
  stop("Could not find src/core/pipeline.R. Please ensure you are running from the correct directory.")
}

# Function to set up argument parser for multi-sample analysis
setup_multi_sample_parser <- function() {
  parser <- ArgumentParser(description = 'Multi-Sample snRNA-seq Analysis Pipeline')
  
  # Input arguments
  parser$add_argument('--rds_inputs',
                      type = 'character',
                      nargs = '+',
                      help = 'Paths to Seurat RDS files to load (multiple files supported).')
  parser$add_argument('--h5_inputs',
                      type = 'character',
                      nargs = '+',
                      help = 'Paths to raw H5 data files (multiple files supported).')
  parser$add_argument('--soupor_cell_doublet_inputs',
                      type = 'character',
                      nargs = '+',
                      help = 'SouporCell clusters output files (multiple files supported).')
  
  # Project settings
  parser$add_argument('--project_name',
                      type = 'character',
                      default = 'MultiSampleProject',
                      help = 'Project name for the multi-sample analysis.')
  parser$add_argument('--working_dir',
                      type = 'character',
                      default = '.',
                      help = 'Working directory where output files will be saved.')
  parser$add_argument('--config_file',
                      type = 'character',
                      default = 'config/settings.yaml',
                      help = 'Path to configuration file.')
  parser$add_argument('--sample_configs',
                      type = 'character',
                      help = 'Path to sample-specific configuration file (YAML).')
  
  # Parallel processing settings
  parser$add_argument('--parallel',
                      type = 'logical',
                      default = TRUE,
                      help = 'Whether to use parallel processing.')
  parser$add_argument('--n_cores',
                      type = 'integer',
                      default = 2,
                      help = 'Number of cores to use for parallel processing (default: 2).')
  
  # Quality control parameters
  parser$add_argument('--min_features',
                      type = 'integer',
                      default = 200,
                      help = 'Minimum number of features per cell for QC filtering.')
  parser$add_argument('--min_counts',
                      type = 'integer',
                      default = 1000,
                      help = 'Minimum number of counts per cell for QC filtering.')
  
  # Processing parameters
  parser$add_argument('--n_variable_features',
                      type = 'integer',
                      default = 2000,
                      help = 'Number of variable features to find.')
  parser$add_argument('--normalization_method',
                      type = 'character',
                      default = 'CLR',
                      choices = c('LogNormalize', 'RC', 'CLR'),
                      help = 'Normalization method to use.')
  parser$add_argument('--scaling_method',
                      type = 'character',
                      default = 'negbinom',
                      choices = c('negbinom', 'linear'),
                      help = 'Scaling method to use.')
  parser$add_argument('--pca_dimensions',
                      type = 'integer',
                      default = 15,
                      help = 'Number of principal components to compute.')
  
  # Clustering parameters
  parser$add_argument('--clustering_resolution',
                      type = 'double',
                      default = 0.5,
                      help = 'Resolution parameter for FindClusters.')
  parser$add_argument('--clustering_algorithm',
                      type = 'character',
                      default = 'leiden',
                      choices = c('louvain', 'multilevel', 'leiden', 'slm'),
                      help = 'Clustering algorithm to use.')
  
  # Analysis options
  parser$add_argument('--find_markers',
                      type = 'logical',
                      default = TRUE,
                      help = 'Whether to find cluster markers.')
  parser$add_argument('--save_intermediate',
                      type = 'logical',
                      default = TRUE,
                      help = 'Whether to save intermediate results.')
  parser$add_argument('--verbose',
                      type = 'logical',
                      default = TRUE,
                      help = 'Whether to print verbose output.')
  
  return(parser)
}

# Main function
main <- function() {
  # Set up argument parser
  parser <- setup_multi_sample_parser()
  args <- parser$parse_args()
  
  # Initialize logger
  init_logger(args$verbose)
  
  # Print welcome message
  cat("=== Multi-Sample snRNA-seq Pipeline ===\n")
  cat("Starting multi-sample analysis...\n")
  cat("Project:", args$project_name, "\n")
  cat("Working directory:", args$working_dir, "\n")
  
  if (!is.null(args$rds_inputs)) {
    cat("RDS inputs:", length(args$rds_inputs), "files\n")
    for (i in seq_along(args$rds_inputs)) {
      cat("  ", i, ":", args$rds_inputs[i], "\n")
    }
  }
  
  if (!is.null(args$h5_inputs)) {
    cat("H5 inputs:", length(args$h5_inputs), "files\n")
    for (i in seq_along(args$h5_inputs)) {
      cat("  ", i, ":", args$h5_inputs[i], "\n")
    }
  }
  
  if (!is.null(args$soupor_cell_doublet_inputs)) {
    cat("SouporCell inputs:", length(args$soupor_cell_doublet_inputs), "files\n")
  }
  
  cat("Parallel processing:", ifelse(args$parallel, "Enabled", "Disabled"), "\n")
  if (args$parallel && !is.null(args$n_cores)) {
    cat("Number of cores:", args$n_cores, "\n")
  }
  cat("\n")
  
  # Run pipeline
  tryCatch({
    result <- run_multi_sample_pipeline(args)
    
    if (!is.null(result)) {
      cat("\n=== Pipeline Completed Successfully ===\n")
      cat("Project:", result$project_name, "\n")
      cat("Total samples:", result$n_samples, "\n")
      cat("Successful samples:", result$n_successful, "\n")
      cat("Failed samples:", result$n_failed, "\n")
      cat("Output directory:", result$project_output_dir, "\n")
      
      if (!is.null(result$summary) && nrow(result$summary) > 0) {
        cat("\nSample Summary:\n")
        print(result$summary)
      }
      
      quit(status = 0)
    } else {
      cat("\nPipeline returned null result\n")
      quit(status = 1)
    }
  }, error = function(e) {
    cat("\nPipeline failed with error:", e$message, "\n")
    quit(status = 1)
  })
}

# Run main function if script is executed directly
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
}

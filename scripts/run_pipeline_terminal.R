#!/usr/bin/env Rscript

# Terminal Interface for snRNA-seq Pipeline
# This script provides a command-line interface for running the pipeline

# Set library path
.libPaths(c("~/.local/lib/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(fs)
  library(future)
  library(promises)
})

# Source pipeline functions
source("src/core/pipeline.R")
source("src/utils/logger.R")
source("src/utils/config.R")

# Set up argument parser for terminal interface
setup_terminal_parser <- function() {
  parser <- ArgumentParser(description = 'snRNA-seq Pipeline Terminal Interface')
  
  # Input options
  parser$add_argument('--rds_input',
                      type = 'character',
                      help = 'Path to a Seurat RDS file to load.')
  parser$add_argument('--h5_input',
                      type = 'character',
                      help = 'Path to a raw H5 data file (e.g., from Cell Ranger).')
  parser$add_argument('--project_name',
                      type = 'character',
                      default = 'SeuratProject',
                      help = 'Project name for the Seurat object.')
  parser$add_argument('--working_dir',
                      type = 'character',
                      default = '.',
                      help = 'Working directory for output files.')
  parser$add_argument('--config_file',
                      type = 'character',
                      default = 'config/settings.yaml',
                      help = 'Path to configuration file.')
  
  # SouporCell options
  parser$add_argument('--soupor_cell_doublet_input',
                      type = 'character',
                      help = 'SouporCell clusters output to identify and remove doublets.')
  
  # Quality control parameters
  parser$add_argument('--min_features',
                      type = 'integer',
                      default = 200,
                      help = 'Minimum number of features per cell for QC filtering.')
  parser$add_argument('--min_counts',
                      type = 'integer',
                      default = 1000,
                      help = 'Minimum number of counts per cell for QC filtering.')
  parser$add_argument('--max_features',
                      type = 'integer',
                      default = 6000,
                      help = 'Maximum number of features per cell for QC filtering.')
  parser$add_argument('--max_counts',
                      type = 'integer',
                      default = 25000,
                      help = 'Maximum number of counts per cell for QC filtering.')
  parser$add_argument('--max_mt_percent',
                      type = 'double',
                      default = 20.0,
                      help = 'Maximum mitochondrial percentage for QC filtering.')
  parser$add_argument('--min_cells',
                      type = 'integer',
                      default = 3,
                      help = 'Minimum number of cells per gene for QC filtering.')
  
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
  parser$add_argument('--scale_factor',
                      type = 'integer',
                      default = 10000,
                      help = 'Scale factor for normalization.')
  
  # Clustering parameters
  parser$add_argument('--clustering_resolution',
                      type = 'double',
                      default = 0.5,
                      help = 'Resolution parameter for clustering.')
  parser$add_argument('--clustering_algorithm',
                      type = 'character',
                      default = 'leiden',
                      choices = c('louvain', 'multilevel', 'leiden', 'slm'),
                      help = 'Clustering algorithm to use.')
  parser$add_argument('--min_cluster_size',
                      type = 'integer',
                      default = 10,
                      help = 'Minimum cluster size.')
  
  # UMAP parameters
  parser$add_argument('--umap_n_neighbors',
                      type = 'integer',
                      default = 30,
                      help = 'UMAP n_neighbors parameter.')
  parser$add_argument('--umap_min_dist',
                      type = 'double',
                      default = 0.3,
                      help = 'UMAP min_dist parameter.')
  
  # Analysis options
  parser$add_argument('--find_markers',
                      action = 'store_true',
                      help = 'Find cluster markers.')
  parser$add_argument('--save_intermediate',
                      action = 'store_true',
                      help = 'Save intermediate results.')
  parser$add_argument('--random_seed',
                      type = 'integer',
                      default = 42,
                      help = 'Random seed for reproducibility.')
  
  # Output options
  parser$add_argument('--verbose',
                      action = 'store_true',
                      help = 'Enable verbose output.')
  parser$add_argument('--output_format',
                      type = 'character',
                      default = 'png',
                      choices = c('png', 'pdf', 'svg'),
                      help = 'Output format for plots.')
  
  return(parser)
}

# Main terminal function
main_terminal <- function() {
  cat("=== snRNA-seq Pipeline Terminal Interface ===\n")
  
  # Set up argument parser
  parser <- setup_terminal_parser()
  args <- parser$parse_args()
  
  # Validate required arguments
  if (is.null(args$rds_input) && is.null(args$h5_input)) {
    cat("Error: Either --rds_input or --h5_input must be specified.\n")
    cat("Use --help for usage information.\n")
    quit(status = 1)
  }
  
  # Initialize logger
  init_logger(args$verbose)
  
  # Convert arguments to list format expected by pipeline
  pipeline_args <- list(
    rds_input = args$rds_input,
    h5_input = args$h5_input,
    project_name = args$project_name,
    working_dir = args$working_dir,
    config_file = args$config_file,
    soupor_cell_doublet_input = args$soupor_cell_doublet_input,
    min_features = args$min_features,
    min_counts = args$min_counts,
    max_features = args$max_features,
    max_counts = args$max_counts,
    max_mt_percent = args$max_mt_percent,
    min_cells = args$min_cells,
    n_variable_features = args$n_variable_features,
    normalization_method = args$normalization_method,
    scaling_method = args$scaling_method,
    pca_dimensions = args$pca_dimensions,
    scale_factor = args$scale_factor,
    clustering_resolution = args$clustering_resolution,
    clustering_algorithm = args$clustering_algorithm,
    min_cluster_size = args$min_cluster_size,
    umap_n_neighbors = args$umap_n_neighbors,
    umap_min_dist = args$umap_min_dist,
    find_markers = args$find_markers,
    save_intermediate = args$save_intermediate,
    random_seed = args$random_seed,
    verbose = args$verbose,
    output_format = args$output_format
  )
  
  # Run pipeline
  cat("Starting pipeline execution...\n")
  tryCatch({
    result <- run_pipeline(pipeline_args)
    if (!is.null(result)) {
      cat("✓ Pipeline completed successfully!\n")
      cat("Output directory:", file.path(args$working_dir, args$project_name), "\n")
      quit(status = 0)
    } else {
      cat("✗ Pipeline failed.\n")
      quit(status = 1)
    }
  }, error = function(e) {
    cat("✗ Pipeline failed with error:", e$message, "\n")
    quit(status = 1)
  })
}

# Run main function if script is executed directly
if (!interactive()) {
  main_terminal()
}

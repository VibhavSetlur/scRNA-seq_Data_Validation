#!/usr/bin/env Rscript

# snRNA-seq Pipeline - Main Script
# This is the main entry point for the pipeline with multi-sample support

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(fs)
  library(Seurat)
  library(future)
  library(future.apply)
  library(parallel)
  library(dplyr)
  library(purrr)
})

# Source utility functions
# Try multiple possible paths for utility files
utility_paths <- c(
  "src/utils/logger.R",
  "../src/utils/logger.R",
  file.path(dirname(getwd()), "src", "utils", "logger.R")
)

logger_found <- FALSE
for (path in utility_paths) {
  if (file.exists(path)) {
    source(path)
    logger_found <- TRUE
    break
  }
}

if (!logger_found) {
  # Create a simple logger if not found
  logger_info <- function(message, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] [%s] %s\n", timestamp, level, message))
  }
  logger_warning <- function(message) logger_info(message, "WARNING")
  logger_error <- function(message) logger_info(message, "ERROR")
  logger_success <- function(message) logger_info(message, "SUCCESS")
}

# Source other utilities
config_paths <- c(
  "src/utils/config.R",
  "../src/utils/config.R",
  file.path(dirname(getwd()), "src", "utils", "config.R")
)

config_found <- FALSE
for (path in config_paths) {
  if (file.exists(path)) {
    source(path)
    config_found <- TRUE
    break
  }
}

if (!config_found) {
  # Create a simple config loader if not found
  load_configuration <- function(config_file) {
    if (file.exists(config_file)) {
      yaml::read_yaml(config_file)
    } else {
      list()
    }
  }
}

validation_paths <- c(
  "src/utils/validation.R",
  "../src/utils/validation.R",
  file.path(dirname(getwd()), "src", "utils", "validation.R")
)

validation_found <- FALSE
for (path in validation_paths) {
  if (file.exists(path)) {
    source(path)
    validation_found <- TRUE
    break
  }
}

if (!validation_found) {
  # Create simple validation functions if not found
  validate_file_exists <- function(file_path) {
    if (!file.exists(file_path)) {
      stop(paste("File not found:", file_path))
    }
  }
}

# Source core modules
# Try multiple possible paths for core modules
core_modules <- c("data_loading.R", "quality_control.R", "processing.R", "clustering.R", "visualization.R")

for (module in core_modules) {
  module_paths <- c(
    file.path("src", "core", module),
    file.path("..", "src", "core", module),
    file.path(dirname(getwd()), "src", "core", module)
  )
  
  module_found <- FALSE
  for (path in module_paths) {
    if (file.exists(path)) {
      source(path)
      module_found <- TRUE
      break
    }
  }
  
  if (!module_found) {
    warning(paste("Could not find module:", module))
  }
}

# Function to set up argument parsing
setup_argument_parser <- function() {
  parser <- ArgumentParser(description = 'snRNA-seq Analysis Pipeline with Multi-Sample Support')
  
  # Input arguments - now supporting multiple files
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
                      default = NULL,
                      help = 'Number of cores to use for parallel processing.')
  
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

# Function to validate arguments for multi-sample analysis
validate_arguments <- function(args) {
  logger_info("Validating arguments for multi-sample analysis...")
  
  # Check that either rds_inputs or h5_inputs is provided
  if (is.null(args$rds_inputs) && is.null(args$h5_inputs)) {
    stop("Either --rds_inputs or --h5_inputs must be provided.")
  }
  
  # Check that not both are provided
  if (!is.null(args$rds_inputs) && !is.null(args$h5_inputs)) {
    stop("Only one of --rds_inputs or --h5_inputs should be provided.")
  }
  
  # Check file existence for all inputs
  if (!is.null(args$rds_inputs)) {
    for (file_path in args$rds_inputs) {
      if (!file.exists(file_path)) {
        stop(paste("RDS input file not found:", file_path))
      }
    }
  }
  
  if (!is.null(args$h5_inputs)) {
    for (file_path in args$h5_inputs) {
      if (!file.exists(file_path)) {
        stop(paste("H5 input file not found:", file_path))
      }
    }
  }
  
  if (!is.null(args$soupor_cell_doublet_inputs)) {
    for (file_path in args$soupor_cell_doublet_inputs) {
      if (!file.exists(file_path)) {
        stop(paste("SouporCell input file not found:", file_path))
      }
    }
  }
  
  # Check configuration file
  if (!file.exists(args$config_file)) {
    stop(paste("Configuration file not found:", args$config_file))
  }
  
  # Check sample configs if provided
  if (!is.null(args$sample_configs) && !file.exists(args$sample_configs)) {
    stop(paste("Sample configuration file not found:", args$sample_configs))
  }
  
  logger_info("Arguments validated successfully.")
}

# Function to create output directory structure for multi-sample analysis
create_output_directory <- function(working_dir, project_name) {
  # Validate inputs
  if (is.null(working_dir) || is.null(project_name)) {
    stop("Working directory and project name cannot be null")
  }
  
  # Ensure working directory exists
  if (!dir.exists(working_dir)) {
    dir.create(working_dir, recursive = TRUE)
    logger_info(paste("Created working directory:", working_dir))
  }
  
  # Create project output directory
  project_output_dir <- file.path(working_dir, paste0(project_name, '_outputs'))
  
  # Validate output directory path
  if (is.null(project_output_dir) || project_output_dir == "" || is.na(project_output_dir)) {
    stop("Invalid output directory path")
  }
  
  if (!dir.exists(project_output_dir)) {
    dir.create(project_output_dir, recursive = TRUE)
    logger_info(paste("Created project output directory:", project_output_dir))
  } else {
    logger_info(paste("Project output directory already exists:", project_output_dir))
  }
  
  # Create subdirectories for different types of outputs
  subdirs <- c("individual_samples", "combined_analysis", "comparisons", "reports")
  for (subdir in subdirs) {
    subdir_path <- file.path(project_output_dir, subdir)
    if (!dir.exists(subdir_path)) {
      dir.create(subdir_path, recursive = TRUE)
      logger_info(paste("Created subdirectory:", subdir_path))
    }
  }
  
  return(project_output_dir)
}

# Function to save session info
save_session_info <- function(output_dir, project_name) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    logger_info(paste("Created output directory for session info:", output_dir))
  }
  
  session_file <- file.path(output_dir, paste0(project_name, '_session_info.txt'))
  
  # Validate file path
  if (is.null(session_file) || session_file == "" || is.na(session_file)) {
    logger_error("Invalid session file path")
    return()
  }
  
  # Ensure the directory exists
  session_dir <- dirname(session_file)
  if (!dir.exists(session_dir)) {
    dir.create(session_dir, recursive = TRUE)
    logger_info(paste("Created session directory:", session_dir))
  }
  
  tryCatch({
    sink(session_file)
    cat("Session Information\n")
    cat("==================\n")
    cat("Date:", Sys.time(), "\n")
    cat("R Version:", R.version.string, "\n")
    cat("Platform:", R.version$platform, "\n")
    cat("\n")
    
    cat("Installed Packages:\n")
    print(installed.packages()[, c("Package", "Version")])
    
    sink()
    
    logger_info(paste("Session info saved to:", session_file))
  }, error = function(e) {
    logger_error(paste("Failed to save session info:", e$message))
    # Ensure sink is closed even if there's an error
    sink()
  })
}

# Function to load sample-specific configurations
load_sample_configurations <- function(sample_configs_file, default_config) {
  if (is.null(sample_configs_file) || !file.exists(sample_configs_file)) {
    logger_info("No sample-specific configuration file provided, using default for all samples")
    return(NULL)
  }
  
  tryCatch({
    sample_configs <- yaml::read_yaml(sample_configs_file)
    logger_info(paste("Loaded sample configurations from:", sample_configs_file))
    return(sample_configs)
  }, error = function(e) {
    logger_warning(paste("Error loading sample configurations:", e$message))
    logger_info("Using default configuration for all samples")
    return(NULL)
  })
}

# Function to get sample-specific configuration
get_sample_config <- function(sample_name, sample_configs, default_config) {
  if (is.null(sample_configs) || is.null(sample_configs[[sample_name]])) {
    return(default_config)
  }
  
  # Merge sample-specific config with default config
  merged_config <- default_config
  sample_config <- sample_configs[[sample_name]]
  
  # Recursively merge configurations
  for (section in names(sample_config)) {
    if (section %in% names(merged_config)) {
      merged_config[[section]] <- modifyList(merged_config[[section]], sample_config[[section]])
    } else {
      merged_config[[section]] <- sample_config[[section]]
    }
  }
  
  return(merged_config)
}

# Function to extract sample name from file path
extract_sample_name <- function(file_path) {
  # Remove directory path and file extension
  base_name <- basename(file_path)
  sample_name <- tools::file_path_sans_ext(base_name)
  
  # Clean up the sample name (remove common suffixes)
  sample_name <- gsub("_filtered_feature_bc_matrix", "", sample_name)
  sample_name <- gsub("_raw_feature_bc_matrix", "", sample_name)
  sample_name <- gsub("_processed", "", sample_name)
  
  return(sample_name)
}

# Function to process a single sample
process_single_sample <- function(sample_info, args, default_config, sample_configs, project_output_dir) {
  sample_name <- sample_info$name
  input_file <- sample_info$file
  input_type <- sample_info$type
  souporcell_file <- sample_info$souporcell_file
  
  logger_info(paste("Processing sample:", sample_name))
  
  tryCatch({
    # Create sample-specific output directory
    sample_output_dir <- file.path(project_output_dir, "individual_samples", sample_name)
    if (!dir.exists(sample_output_dir)) {
      dir.create(sample_output_dir, recursive = TRUE)
    }
    
    # Get sample-specific configuration
    sample_config <- get_sample_config(sample_name, sample_configs, default_config)
    
    # Create sample-specific arguments
    sample_args <- args
    sample_args$project_name <- sample_name
    sample_args$working_dir <- sample_output_dir
    
    # Set input file based on type
    if (input_type == "rds") {
      sample_args$rds_input <- input_file
      sample_args$h5_input <- NULL
    } else {
      sample_args$h5_input <- input_file
      sample_args$rds_input <- NULL
    }
    
    # Set SouporCell file if available
    if (!is.null(souporcell_file)) {
      sample_args$soupor_cell_doublet_input <- souporcell_file
    }
    
    # Run pipeline for this sample
    result <- run_single_sample_pipeline(sample_args, sample_config, sample_output_dir)
    
    # Add sample metadata
    result$sample_name <- sample_name
    result$input_file <- input_file
    result$input_type <- input_type
    
    logger_success(paste("Sample", sample_name, "processed successfully"))
    return(result)
    
  }, error = function(e) {
    logger_error(paste("Error processing sample", sample_name, ":", e$message))
    return(list(
      sample_name = sample_name,
      error = e$message,
      status = "failed"
    ))
  })
}

# Function to run pipeline for a single sample
run_single_sample_pipeline <- function(args, config, output_dir) {
  logger_info(paste("Running pipeline for sample:", args$project_name))
  
  tryCatch({
    # Load data
    seurat_object <- load_data(args, config)
    
    if (is.null(seurat_object) || ncol(seurat_object) == 0) {
      stop("No valid data loaded for this sample")
    }
    
    # Quality control
    seurat_object <- perform_quality_control(seurat_object, args, config, output_dir)
    
    if (is.null(seurat_object) || ncol(seurat_object) == 0) {
      stop("No cells remaining after quality control for this sample")
    }
    
    # Data processing
    seurat_object <- process_data(seurat_object, args, config, output_dir)
    
    if (is.null(seurat_object)) {
      stop("Data processing failed for this sample")
    }
    
    # Clustering and dimensionality reduction
    seurat_object <- perform_clustering(seurat_object, args, config, output_dir)
    
    if (is.null(seurat_object)) {
      stop("Clustering failed for this sample")
    }
    
    # Generate visualizations
    generate_visualizations(seurat_object, args, config, output_dir)
    
    # Save final object
    final_file <- file.path(output_dir, paste0(args$project_name, '_processed.rds'))
    saveRDS(seurat_object, file = final_file)
    
    # Collect results
    pipeline_results <- list(
      project_name = args$project_name,
      output_dir = output_dir,
      n_cells = ncol(seurat_object),
      n_genes = nrow(seurat_object),
      n_clusters = length(unique(Seurat::Idents(seurat_object))),
      seurat_object = seurat_object,
      final_file = final_file,
      status = "completed"
    )
    
    return(pipeline_results)
    
  }, error = function(e) {
    logger_error(paste("Pipeline failed for sample", args$project_name, ":", e$message))
    return(list(
      project_name = args$project_name,
      error = e$message,
      status = "failed"
    ))
  })
}

# Function to set up parallel processing
setup_parallel_processing <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- min(parallel::detectCores(), 2)  # Conservative approach, max 2
  }
  
  # Ensure n_cores doesn't exceed available cores and respect system limits
  available_cores <- parallel::detectCores()
  # Use conservative approach - don't exceed 75% of available cores
  max_safe_cores <- max(1, floor(available_cores * 0.75))
  n_cores <- min(n_cores, max_safe_cores, 2)  # Cap at 2 cores for safety
  
  logger_info(paste("Setting up parallel processing with", n_cores, "cores"))
  
  # Set up future plan with error handling
  if (n_cores > 1) {
    tryCatch({
      plan(multisession, workers = n_cores)
      logger_success(paste("Parallel processing enabled with", n_cores, "workers"))
    }, error = function(e) {
      logger_warning(paste("Parallel processing failed, falling back to sequential:", e$message))
      plan(sequential)
      n_cores <- 1
    })
  } else {
    plan(sequential)
    logger_info("Using sequential processing")
  }
  
  return(n_cores)
}

# Function to create sample information list
create_sample_info_list <- function(args) {
  sample_info_list <- list()
  
  # Process RDS inputs
  if (!is.null(args$rds_inputs)) {
    for (i in seq_along(args$rds_inputs)) {
      sample_name <- extract_sample_name(args$rds_inputs[i])
      sample_info <- list(
        name = sample_name,
        file = args$rds_inputs[i],
        type = "rds",
        souporcell_file = if (!is.null(args$soupor_cell_doublet_inputs) && i <= length(args$soupor_cell_doublet_inputs)) 
                           args$soupor_cell_doublet_inputs[i] else NULL
      )
      sample_info_list[[sample_name]] <- sample_info
    }
  }
  
  # Process H5 inputs
  if (!is.null(args$h5_inputs)) {
    for (i in seq_along(args$h5_inputs)) {
      sample_name <- extract_sample_name(args$h5_inputs[i])
      sample_info <- list(
        name = sample_name,
        file = args$h5_inputs[i],
        type = "h5",
        souporcell_file = if (!is.null(args$soupor_cell_doublet_inputs) && i <= length(args$soupor_cell_doublet_inputs)) 
                           args$soupor_cell_doublet_inputs[i] else NULL
      )
      sample_info_list[[sample_name]] <- sample_info
    }
  }
  
  return(sample_info_list)
}

# Function to generate multi-sample summary report
generate_multi_sample_summary <- function(results, project_output_dir, project_name) {
  logger_info("Generating multi-sample summary report...")
  
  # Create summary data
  summary_data <- data.frame(
    Sample_Name = character(),
    Status = character(),
    Cells = integer(),
    Genes = integer(),
    Clusters = integer(),
    stringsAsFactors = FALSE
  )
  
  successful_samples <- list()
  
  for (sample_name in names(results)) {
    result <- results[[sample_name]]
    
    if (result$status == "completed") {
      summary_data <- rbind(summary_data, data.frame(
        Sample_Name = sample_name,
        Status = "Completed",
        Cells = result$n_cells,
        Genes = result$n_genes,
        Clusters = result$n_clusters,
        stringsAsFactors = FALSE
      ))
      successful_samples[[sample_name]] <- result
    } else {
      summary_data <- rbind(summary_data, data.frame(
        Sample_Name = sample_name,
        Status = "Failed",
        Cells = 0,
        Genes = 0,
        Clusters = 0,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save summary
  summary_file <- file.path(project_output_dir, paste0(project_name, '_multi_sample_summary.csv'))
  write.csv(summary_data, file = summary_file, row.names = FALSE)
  
  # Generate summary plots
  if (length(successful_samples) > 0) {
    generate_multi_sample_plots(successful_samples, project_output_dir, project_name)
  }
  
  logger_success(paste("Multi-sample summary saved to:", summary_file))
  return(summary_data)
}

# Function to generate multi-sample comparison plots
generate_multi_sample_plots <- function(successful_samples, project_output_dir, project_name) {
  logger_info("Generating multi-sample comparison plots...")
  
  # Create comparison directory
  comparison_dir <- file.path(project_output_dir, "comparisons")
  if (!dir.exists(comparison_dir)) {
    dir.create(comparison_dir, recursive = TRUE)
  }
  
  # Sample comparison plot
  sample_comparison_data <- data.frame(
    Sample = names(successful_samples),
    Cells = sapply(successful_samples, function(x) x$n_cells),
    Genes = sapply(successful_samples, function(x) x$n_genes),
    Clusters = sapply(successful_samples, function(x) x$n_clusters)
  )
  
  # Create comparison plots
  p1 <- ggplot2::ggplot(sample_comparison_data, ggplot2::aes(x = Sample, y = Cells, fill = Sample)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0(project_name, " - Cell Counts by Sample")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot2::ggplot(sample_comparison_data, ggplot2::aes(x = Sample, y = Genes, fill = Sample)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0(project_name, " - Gene Counts by Sample")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot2::ggplot(sample_comparison_data, ggplot2::aes(x = Sample, y = Clusters, fill = Sample)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0(project_name, " - Cluster Counts by Sample")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Save plots
  ggplot2::ggsave(
    filename = file.path(comparison_dir, paste0(project_name, '_cell_counts_comparison.png')),
    plot = p1,
    width = 10,
    height = 6,
    units = 'in',
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(comparison_dir, paste0(project_name, '_gene_counts_comparison.png')),
    plot = p2,
    width = 10,
    height = 6,
    units = 'in',
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(comparison_dir, paste0(project_name, '_cluster_counts_comparison.png')),
    plot = p3,
    width = 10,
    height = 6,
    units = 'in',
    dpi = 300
  )
  
  logger_success("Multi-sample comparison plots generated")
}

# Main pipeline function for multi-sample analysis
run_multi_sample_pipeline <- function(args) {
  logger_info("=== Starting Multi-Sample snRNA-seq Pipeline ===")
  
  tryCatch({
    # Validate arguments
    validate_arguments(args)
    
    # Load configuration
    if (!is.null(args$config_file) && file.exists(args$config_file)) {
      default_config <- load_configuration(args$config_file)
    } else {
      logger_info("No configuration file provided, using default settings")
      default_config <- create_default_config()
    }
    
    # Load sample-specific configurations
    sample_configs <- load_sample_configurations(args$sample_configs, default_config)
    
    # Create output directory
    project_output_dir <- create_output_directory(args$working_dir, args$project_name)
    
    # Save session info
    save_session_info(project_output_dir, args$project_name)
    
    # Set up parallel processing
    n_cores <- setup_parallel_processing(args$n_cores)
    
    # Create sample information list
    sample_info_list <- create_sample_info_list(args)
    
    logger_info(paste("Processing", length(sample_info_list), "samples"))
    
    # Process samples in parallel
    if (args$parallel && n_cores > 1) {
      logger_info("Processing samples in parallel...")
      results <- future_lapply(
        sample_info_list,
        function(sample_info) {
          process_single_sample(sample_info, args, default_config, sample_configs, project_output_dir)
        }
      )
    } else {
      logger_info("Processing samples sequentially...")
      results <- lapply(
        sample_info_list,
        function(sample_info) {
          process_single_sample(sample_info, args, default_config, sample_configs, project_output_dir)
        }
      )
    }
    
    # Generate multi-sample summary
    summary_data <- generate_multi_sample_summary(results, project_output_dir, args$project_name)
    
    # Create final results object
    final_results <- list(
      project_name = args$project_name,
      project_output_dir = project_output_dir,
      samples = results,
      summary = summary_data,
      n_samples = length(sample_info_list),
      n_successful = sum(sapply(results, function(x) x$status == "completed")),
      n_failed = sum(sapply(results, function(x) x$status == "failed")),
      parallel_cores = n_cores
    )
    
    logger_success("=== Multi-Sample Pipeline completed successfully ===")
    return(final_results)
    
  }, error = function(e) {
    logger_error(paste("Multi-sample pipeline failed with error:", e$message))
    stop(paste("Multi-sample pipeline execution failed:", e$message))
  })
}

# Main execution
main <- function() {
  # Set up argument parser
  parser <- setup_argument_parser()
  args <- parser$parse_args()
  
  # Initialize logger
  init_logger(args$verbose)
  
  # Run pipeline
  tryCatch({
    result <- run_multi_sample_pipeline(args)
    if (!is.null(result)) {
      quit(status = 0)
    } else {
      quit(status = 1)
    }
  }, error = function(e) {
    logger_error(paste("Pipeline failed with error:", e$message))
    quit(status = 1)
  })
}

# Run main function if script is executed directly
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
}

# Validation Utility for snRNA-seq Pipeline
# Provides input validation functions

# Function to validate file paths
validate_file_paths <- function(args) {
  logger_info("Validating file paths...")
  
  # Check input files
  if (!is.null(args$rds_input)) {
    if (!file.exists(args$rds_input)) {
      stop(paste("RDS input file not found:", args$rds_input))
    }
  }
  
  if (!is.null(args$h5_input)) {
    if (!file.exists(args$h5_input)) {
      stop(paste("H5 input file not found:", args$h5_input))
    }
  }
  
  if (!is.null(args$soupor_cell_doublet_input)) {
    if (!file.exists(args$soupor_cell_doublet_input)) {
      stop(paste("SouporCell input file not found:", args$soupor_cell_doublet_input))
    }
  }
  
  # Check configuration file
  if (!file.exists(args$config_file)) {
    stop(paste("Configuration file not found:", args$config_file))
  }
  
  logger_success("File path validation passed")
}

# Function to validate parameters
validate_parameters <- function(args) {
  logger_info("Validating parameters...")
  
  # Validate numeric parameters
  numeric_params <- list(
    min_features = args$min_features,
    min_counts = args$min_counts,
    n_variable_features = args$n_variable_features,
    pca_dimensions = args$pca_dimensions,
    clustering_resolution = args$clustering_resolution
  )
  
  for (param_name in names(numeric_params)) {
    value <- numeric_params[[param_name]]
    if (!is.numeric(value) || value <= 0) {
      stop(paste("Parameter", param_name, "must be a positive number"))
    }
  }
  
  # Validate clustering resolution
  if (args$clustering_resolution <= 0 || args$clustering_resolution > 2) {
    stop("Clustering resolution must be between 0 and 2")
  }
  
  # Validate clustering algorithm
  valid_algorithms <- c("louvain", "multilevel", "leiden", "slm")
  if (!args$clustering_algorithm %in% valid_algorithms) {
    stop(paste("Invalid clustering algorithm. Must be one of:", paste(valid_algorithms, collapse = ", ")))
  }
  
  # Validate normalization method
  valid_norm_methods <- c("LogNormalize", "RC", "CLR")
  if (!args$normalization_method %in% valid_norm_methods) {
    stop(paste("Invalid normalization method. Must be one of:", paste(valid_norm_methods, collapse = ", ")))
  }
  
  # Validate scaling method
  valid_scaling_methods <- c("negbinom", "linear")
  if (!args$scaling_method %in% valid_scaling_methods) {
    stop(paste("Invalid scaling method. Must be one of:", paste(valid_scaling_methods, collapse = ", ")))
  }
  
  logger_success("Parameter validation passed")
}

# Function to validate system requirements
validate_system_requirements <- function() {
  logger_info("Validating system requirements...")
  
  # Check R version
  r_version <- R.version.string
  major_version <- as.numeric(R.version$major)
  minor_version <- as.numeric(R.version$minor)
  
  if (major_version < 4 || (major_version == 4 && minor_version < 3)) {
    logger_warning("R version 4.3.0 or higher is recommended")
  }
  
  # Check available memory
  if (requireNamespace("pryr", quietly = TRUE)) {
    mem_available <- pryr::mem_used()
    if (mem_available < 4 * 1024^3) {  # 4GB
      logger_warning("Less than 4GB of memory available. Large datasets may cause issues.")
    }
  }
  
  # Check disk space
  disk_space <- file.size(".") / (1024^3)  # GB
  if (disk_space < 1) {
    logger_warning("Less than 1GB of disk space available in current directory")
  }
  
  logger_success("System requirements validation completed")
}

# Function to validate Seurat object
validate_seurat_object <- function(seurat_object) {
  logger_info("Validating Seurat object...")
  
  # Check if it's a Seurat object
  if (!inherits(seurat_object, "Seurat")) {
    stop("Object is not a Seurat object")
  }
  
  # Check dimensions
  n_cells <- ncol(seurat_object)
  n_genes <- nrow(seurat_object)
  
  if (n_cells == 0) {
    stop("Seurat object contains no cells")
  }
  
  if (n_genes == 0) {
    stop("Seurat object contains no genes")
  }
  
  # Check for required metadata
  required_meta <- c("nFeature_RNA", "nCount_RNA")
  missing_meta <- setdiff(required_meta, colnames(seurat_object@meta.data))
  
  if (length(missing_meta) > 0) {
    logger_warning(paste("Missing metadata columns:", paste(missing_meta, collapse = ", ")))
  }
  
  logger_success(paste("Seurat object validation passed. Dimensions:", n_genes, "genes x", n_cells, "cells"))
  
  return(TRUE)
}

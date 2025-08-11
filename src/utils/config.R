# Configuration utility for snRNA-seq Pipeline
# Handles loading and validation of configuration files

# Function to load configuration from YAML file
load_configuration <- function(config_file) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("yaml package is required to load configuration")
  }
  
  if (!file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file))
  }
  
  tryCatch({
    config <- yaml::read_yaml(config_file)
    logger_info(paste("Configuration loaded from:", config_file))
    return(config)
  }, error = function(e) {
    stop(paste("Error loading configuration:", e$message))
  })
}

# Function to validate configuration
validate_configuration <- function(config) {
  logger_info("Validating configuration...")
  
  # Required sections
  required_sections <- c("project", "qc", "processing", "clustering", "visualization")
  missing_sections <- setdiff(required_sections, names(config))
  
  if (length(missing_sections) > 0) {
    stop(paste("Missing required configuration sections:", 
               paste(missing_sections, collapse = ", ")))
  }
  
  # Validate project settings
  validate_project_config(config$project)
  
  # Validate QC settings
  validate_qc_config(config$qc)
  
  # Validate processing settings
  validate_processing_config(config$processing)
  
  # Validate clustering settings
  validate_clustering_config(config$clustering)
  
  # Validate visualization settings
  validate_visualization_config(config$visualization)
  
  logger_info("Configuration validation passed")
  return(TRUE)
}

# Function to validate project configuration
validate_project_config <- function(project_config) {
  if (!is.list(project_config)) {
    stop("Project configuration must be a list")
  }
  
  # Check required fields
  if (is.null(project_config$default_name)) {
    stop("Project configuration missing 'default_name'")
  }
  
  if (is.null(project_config$output_format)) {
    stop("Project configuration missing 'output_format'")
  }
  
  # Validate output format
  valid_formats <- c("png", "pdf", "svg")
  if (!project_config$output_format %in% valid_formats) {
    stop(paste("Invalid output format. Must be one of:", 
               paste(valid_formats, collapse = ", ")))
  }
  
  # Validate DPI
  if (!is.null(project_config$dpi)) {
    if (!is.numeric(project_config$dpi) || project_config$dpi <= 0) {
      stop("DPI must be a positive number")
    }
  }
}

# Function to validate QC configuration
validate_qc_config <- function(qc_config) {
  if (!is.list(qc_config)) {
    stop("QC configuration must be a list")
  }
  
  # Check required fields
  required_fields <- c("min_features", "min_counts", "max_features", "max_counts", "max_mt_percent")
  missing_fields <- setdiff(required_fields, names(qc_config))
  
  if (length(missing_fields) > 0) {
    stop(paste("QC configuration missing required fields:", 
               paste(missing_fields, collapse = ", ")))
  }
  
  # Validate numeric fields
  numeric_fields <- c("min_features", "min_counts", "max_features", "max_counts", "max_mt_percent")
  for (field in numeric_fields) {
    if (!is.numeric(qc_config[[field]]) || qc_config[[field]] < 0) {
      stop(paste("QC configuration field", field, "must be a non-negative number"))
    }
  }
  
  # Validate logical relationships
  if (qc_config$min_features > qc_config$max_features) {
    stop("min_features cannot be greater than max_features")
  }
  
  if (qc_config$min_counts > qc_config$max_counts) {
    stop("min_counts cannot be greater than max_counts")
  }
  
  if (qc_config$max_mt_percent > 100) {
    stop("max_mt_percent cannot be greater than 100")
  }
}

# Function to validate processing configuration
validate_processing_config <- function(processing_config) {
  if (!is.list(processing_config)) {
    stop("Processing configuration must be a list")
  }
  
  # Check required fields
  required_fields <- c("normalization_method", "n_variable_features", "scaling_method", "pca_dimensions")
  missing_fields <- setdiff(required_fields, names(processing_config))
  
  if (length(missing_fields) > 0) {
    stop(paste("Processing configuration missing required fields:", 
               paste(missing_fields, collapse = ", ")))
  }
  
  # Validate normalization method
  valid_norm_methods <- c("LogNormalize", "RC", "CLR")
  if (!processing_config$normalization_method %in% valid_norm_methods) {
    stop(paste("Invalid normalization method. Must be one of:", 
               paste(valid_norm_methods, collapse = ", ")))
  }
  
  # Validate scaling method
  valid_scaling_methods <- c("negbinom", "linear")
  if (!processing_config$scaling_method %in% valid_scaling_methods) {
    stop(paste("Invalid scaling method. Must be one of:", 
               paste(valid_scaling_methods, collapse = ", ")))
  }
  
  # Validate numeric fields
  if (!is.numeric(processing_config$n_variable_features) || 
      processing_config$n_variable_features <= 0) {
    stop("n_variable_features must be a positive number")
  }
  
  if (!is.numeric(processing_config$pca_dimensions) || 
      processing_config$pca_dimensions <= 0) {
    stop("pca_dimensions must be a positive number")
  }
}

# Function to validate clustering configuration
validate_clustering_config <- function(clustering_config) {
  if (!is.list(clustering_config)) {
    stop("Clustering configuration must be a list")
  }
  
  # Check required fields
  required_fields <- c("resolution", "algorithm", "min_cluster_size")
  missing_fields <- setdiff(required_fields, names(clustering_config))
  
  if (length(missing_fields) > 0) {
    stop(paste("Clustering configuration missing required fields:", 
               paste(missing_fields, collapse = ", ")))
  }
  
  # Validate algorithm
  valid_algorithms <- c("louvain", "multilevel", "leiden", "slm")
  if (!clustering_config$algorithm %in% valid_algorithms) {
    stop(paste("Invalid clustering algorithm. Must be one of:", 
               paste(valid_algorithms, collapse = ", ")))
  }
  
  # Validate numeric fields
  if (!is.numeric(clustering_config$resolution) || 
      clustering_config$resolution <= 0) {
    stop("Clustering resolution must be a positive number")
  }
  
  if (!is.numeric(clustering_config$min_cluster_size) || 
      clustering_config$min_cluster_size <= 0) {
    stop("min_cluster_size must be a positive number")
  }
}

# Function to validate visualization configuration
validate_visualization_config <- function(viz_config) {
  if (!is.list(viz_config)) {
    stop("Visualization configuration must be a list")
  }
  
  # Check required fields
  required_fields <- c("theme", "color_palette", "point_size", "label_size")
  missing_fields <- setdiff(required_fields, names(viz_config))
  
  if (length(missing_fields) > 0) {
    stop(paste("Visualization configuration missing required fields:", 
               paste(missing_fields, collapse = ", ")))
  }
  
  # Validate theme
  valid_themes <- c("classic", "minimal", "bw")
  if (!viz_config$theme %in% valid_themes) {
    stop(paste("Invalid theme. Must be one of:", 
               paste(valid_themes, collapse = ", ")))
  }
  
  # Validate color palette
  valid_palettes <- c("viridis", "rainbow", "brewer")
  if (!viz_config$color_palette %in% valid_palettes) {
    stop(paste("Invalid color palette. Must be one of:", 
               paste(valid_palettes, collapse = ", ")))
  }
  
  # Validate numeric fields
  numeric_fields <- c("point_size", "label_size", "title_size", "subtitle_size")
  for (field in numeric_fields) {
    if (!is.null(viz_config[[field]])) {
      if (!is.numeric(viz_config[[field]]) || viz_config[[field]] <= 0) {
        stop(paste("Visualization configuration field", field, "must be a positive number"))
      }
    }
  }
}

# Function to merge configuration with command line arguments
merge_config_with_args <- function(config, args) {
  logger_info("Merging configuration with command line arguments...")
  
  # Create merged configuration
  merged_config <- config
  
  # Override with command line arguments if provided
  if (!is.null(args$min_features)) {
    merged_config$qc$min_features <- args$min_features
  }
  
  if (!is.null(args$min_counts)) {
    merged_config$qc$min_counts <- args$min_counts
  }
  
  if (!is.null(args$n_variable_features)) {
    merged_config$processing$n_variable_features <- args$n_variable_features
  }
  
  if (!is.null(args$normalization_method)) {
    merged_config$processing$normalization_method <- args$normalization_method
  }
  
  if (!is.null(args$scaling_method)) {
    merged_config$processing$scaling_method <- args$scaling_method
  }
  
  if (!is.null(args$pca_dimensions)) {
    merged_config$processing$pca_dimensions <- args$pca_dimensions
  }
  
  if (!is.null(args$clustering_resolution)) {
    merged_config$clustering$resolution <- args$clustering_resolution
  }
  
  if (!is.null(args$clustering_algorithm)) {
    merged_config$clustering$algorithm <- args$clustering_algorithm
  }
  
  logger_info("Configuration merged successfully")
  return(merged_config)
}

# Function to save configuration
save_configuration <- function(config, output_file) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("yaml package is required to save configuration")
  }
  
  tryCatch({
    yaml::write_yaml(config, output_file)
    logger_info(paste("Configuration saved to:", output_file))
  }, error = function(e) {
    stop(paste("Error saving configuration:", e$message))
  })
}

# Function to create default configuration
create_default_config <- function() {
  config <- list(
    project = list(
      default_name = "SeuratProject",
      output_format = "png",
      dpi = 300,
      save_intermediate = TRUE
    ),
    qc = list(
      min_features = 200,
      min_counts = 1000,
      max_features = 6000,
      max_counts = 25000,
      max_mt_percent = 20,
      min_cells = 3
    ),
    processing = list(
      normalization_method = "CLR",
      n_variable_features = 2000,
      scaling_method = "negbinom",
      pca_dimensions = 15,
      scale_factor = 10000
    ),
    clustering = list(
      resolution = 0.5,
      algorithm = "leiden",
      min_cluster_size = 10,
      random_seed = 42
    ),
    umap = list(
      n_neighbors = 30,
      min_dist = 0.3,
      spread = 1.0,
      random_seed = 42
    ),
    differential_expression = list(
      find_markers = TRUE,
      test_use = "DESeq2",
      logfc_threshold = 0.25,
      min_pct = 0.1,
      return_thresh = 0.05,
      only_pos = FALSE
    ),
    visualization = list(
      theme = "classic",
      color_palette = "viridis",
      point_size = 0.7,
      label_size = 3,
      title_size = 14,
      subtitle_size = 10
    ),
    performance = list(
      parallel_processing = TRUE,
      n_cores = "auto",
      memory_limit = "8GB",
      chunk_size = 1000
    ),
    paths = list(
      working_directory = ".",
      output_directory = "outputs",
      log_directory = "logs",
      temp_directory = "temp"
    ),
    advanced = list(
      verbose = TRUE,
      save_logs = TRUE,
      check_memory = TRUE,
      validate_inputs = TRUE,
      backup_results = TRUE
    )
  )
  
  return(config)
}

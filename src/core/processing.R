# Data Processing Module for snRNA-seq Pipeline
# Handles normalization, variable feature identification, and scaling

# Function to process data
process_data <- function(seurat_object, args, config, output_dir) {
  logger_info("Starting data processing...")
  
  # Normalize data
  seurat_object <- normalize_data(seurat_object, args, config)
  
  # Find variable features
  seurat_object <- find_variable_features(seurat_object, args, config, output_dir)
  
  # Scale data
  seurat_object <- scale_data(seurat_object, args, config)
  
  logger_success("Data processing completed")
  
  return(seurat_object)
}

# Function to normalize data
normalize_data <- function(seurat_object, args, config) {
  logger_info(paste("Normalizing data using", args$normalization_method, "method..."))
  
  tryCatch({
    # Set scale factor based on normalization method
    scale_factor <- ifelse(args$normalization_method == "CLR", 10000, NULL)
    
    seurat_object <- Seurat::NormalizeData(
      object = seurat_object,
      normalization.method = args$normalization_method,
      scale.factor = scale_factor,
      margin = 2  # cells
    )
    
    logger_success("Data normalization completed")
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error during normalization:", e$message))
    stop(e)
  })
}

# Function to find variable features
find_variable_features <- function(seurat_object, args, config, output_dir) {
  logger_info(paste("Finding", args$n_variable_features, "variable features..."))
  
  tryCatch({
    # Check if there are enough features
    n_features <- nrow(seurat_object)
    n_variable <- min(args$n_variable_features, n_features)
    
    if (n_features < args$n_variable_features) {
      logger_warning(paste("Requested", args$n_variable_features, "variable features, but only", n_features, "features available"))
    }
    
    seurat_object <- Seurat::FindVariableFeatures(
      object = seurat_object,
      selection.method = 'vst',
      nfeatures = n_variable
    )
    
    # Generate variable features plot
    generate_variable_features_plot(seurat_object, args, config, output_dir)
    
    logger_success(paste("Found", length(Seurat::VariableFeatures(seurat_object)), "variable features"))
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error finding variable features:", e$message))
    stop(e)
  })
}

# Function to scale data
scale_data <- function(seurat_object, args, config) {
  logger_info(paste("Scaling data using", args$scaling_method, "model..."))
  
  tryCatch({
    seurat_object <- Seurat::ScaleData(
      object = seurat_object,
      features = Seurat::VariableFeatures(seurat_object),
      do.scale = TRUE,
      do.center = TRUE,
      model.use = args$scaling_method,
      assay = 'RNA'
    )
    
    logger_success("Data scaling completed")
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error during scaling:", e$message))
    stop(e)
  })
}

# Function to generate variable features plot
generate_variable_features_plot <- function(seurat_object, args, config, output_dir) {
  logger_info("Generating variable features plot...")
  
  tryCatch({
    # Get top variable features
    variable_features <- Seurat::VariableFeatures(seurat_object)
    n_features_to_label <- min(20, length(variable_features))
    top_features <- head(variable_features, n_features_to_label)
    
    # Create plot
    p <- Seurat::VariableFeaturePlot(seurat_object) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = paste0(args$project_name, ' - Variable Features (Top ', n_features_to_label, ' Labeled)'),
        x = 'Average Expression',
        y = 'Variance (vst)'
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14))
    
    # Add labels
    if (length(top_features) > 0) {
      p <- Seurat::LabelPoints(
        plot = p,
        points = top_features,
        repel = TRUE,
        size = 3
      )
    }
    
    # Save plot
    plot_file <- file.path(output_dir, paste0(args$project_name, '_Variable_Features.png'))
    ggplot2::ggsave(
      filename = plot_file,
      plot = p,
      width = 10,
      height = 8,
      units = 'in',
      dpi = config$project$dpi
    )
    
    logger_info(paste("Variable features plot saved to:", plot_file))
    
  }, error = function(e) {
    logger_warning(paste("Error generating variable features plot:", e$message))
  })
}

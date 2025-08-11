# Clustering Module for snRNA-seq Pipeline
# Handles PCA, clustering, and UMAP dimensionality reduction

# Function to perform clustering
perform_clustering <- function(seurat_object, args, config, output_dir) {
  logger_info("Starting clustering analysis...")
  
  # Run PCA
  seurat_object <- run_pca(seurat_object, args, config, output_dir)
  
  # Find neighbors
  seurat_object <- find_neighbors(seurat_object, args, config)
  
  # Find clusters
  seurat_object <- find_clusters(seurat_object, args, config)
  
  # Run UMAP
  seurat_object <- run_umap(seurat_object, args, config, output_dir)
  
  # Find markers if requested
  if (args$find_markers) {
    find_markers(seurat_object, args, config, output_dir)
  }
  
  logger_success("Clustering analysis completed")
  
  return(seurat_object)
}

# Function to run PCA
run_pca <- function(seurat_object, args, config, output_dir) {
  logger_info(paste("Running PCA with", args$pca_dimensions, "dimensions..."))
  
  tryCatch({
    # Check if there are enough cells
    if (ncol(seurat_object) <= 1) {
      stop("Not enough cells to perform PCA")
    }
    
    # Check if there are enough variable features
    n_variable <- length(Seurat::VariableFeatures(seurat_object))
    n_pca <- min(args$pca_dimensions, n_variable, ncol(seurat_object) - 1)
    
    if (n_pca < args$pca_dimensions) {
      logger_warning(paste("Requested", args$pca_dimensions, "PCA dimensions, but only", n_pca, "available"))
    }
    
    seurat_object <- Seurat::RunPCA(
      object = seurat_object,
      npcs = n_pca,
      rev.pca = TRUE,
      weight.by.var = TRUE,
      verbose = FALSE
    )
    
    # Generate PCA plots
    generate_pca_plots(seurat_object, args, config, output_dir)
    
    logger_success(paste("PCA completed with", n_pca, "dimensions"))
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error during PCA:", e$message))
    stop(e)
  })
}

# Function to find neighbors
find_neighbors <- function(seurat_object, args, config) {
  logger_info("Finding neighbors...")
  
  tryCatch({
    # Get valid PCA dimensions
    pca_dims <- 1:min(args$pca_dimensions, ncol(seurat_object@reductions$pca@feature.loadings))
    
    seurat_object <- Seurat::FindNeighbors(
      object = seurat_object,
      dims = pca_dims,
      nn.method = 'annoy',
      verbose = FALSE
    )
    
    logger_success("Neighbors found")
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error finding neighbors:", e$message))
    stop(e)
  })
}

# Function to find clusters
find_clusters <- function(seurat_object, args, config) {
  logger_info(paste("Finding clusters using", args$clustering_algorithm, "algorithm..."))
  
  tryCatch({
    # Map algorithm names to numbers
    algorithm_map <- c('louvain' = 1, 'multilevel' = 2, 'slm' = 3, 'leiden' = 4)
    algorithm_num <- algorithm_map[args$clustering_algorithm]
    
    if (is.na(algorithm_num)) {
      logger_warning("Invalid clustering algorithm, using leiden")
      algorithm_num <- 4
    }
    
    seurat_object <- Seurat::FindClusters(
      object = seurat_object,
      resolution = args$clustering_resolution,
      algorithm = algorithm_num,
      verbose = FALSE
    )
    
    logger_success(paste("Found", length(unique(Seurat::Idents(seurat_object))), "clusters"))
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error finding clusters:", e$message))
    stop(e)
  })
}

# Function to run UMAP
run_umap <- function(seurat_object, args, config, output_dir) {
  logger_info("Running UMAP...")
  
  tryCatch({
    # Get valid PCA dimensions
    pca_dims <- 1:min(args$pca_dimensions, ncol(seurat_object@reductions$pca@feature.loadings))
    
    seurat_object <- Seurat::RunUMAP(
      object = seurat_object,
      dims = pca_dims,
      verbose = FALSE
    )
    
    # Generate UMAP plot
    generate_umap_plot(seurat_object, args, config, output_dir)
    
    logger_success("UMAP completed")
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error during UMAP:", e$message))
    stop(e)
  })
}

# Function to find markers
find_markers <- function(seurat_object, args, config, output_dir) {
  logger_info("Finding cluster markers...")
  
  tryCatch({
    # Check if there are multiple clusters
    n_clusters <- length(unique(Seurat::Idents(seurat_object)))
    
    if (n_clusters <= 1) {
      logger_warning("Only one cluster found, skipping marker analysis")
      return(seurat_object)
    }
    
    # Find markers
    markers <- Seurat::FindAllMarkers(
      object = seurat_object,
      logfc.threshold = config$differential_expression$logfc_threshold,
      min.pct = config$differential_expression$min_pct,
      test.use = config$differential_expression$test_use,
      return.thresh = config$differential_expression$return_thresh,
      verbose = FALSE
    )
    
    # Save markers
    markers_file <- file.path(output_dir, paste0(args$project_name, '_markers.csv'))
    write.csv(markers, file = markers_file, row.names = FALSE)
    
    logger_success(paste("Markers saved to:", markers_file))
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error finding markers:", e$message))
    return(seurat_object)
  })
}

# Function to generate PCA plots
generate_pca_plots <- function(seurat_object, args, config, output_dir) {
  logger_info("Generating PCA plots...")
  
  tryCatch({
    # PCA DimPlot
    p1 <- Seurat::DimPlot(
      object = seurat_object,
      reduction = 'pca'
    ) +
      ggplot2::labs(
        title = paste0(args$project_name, ' - PCA Plot'),
        subtitle = paste0('Computed Dimensions: ', args$pca_dimensions)
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      )
    
    # Elbow plot
    p2 <- Seurat::ElbowPlot(
      object = seurat_object,
      reduction = 'pca',
      ndims = args$pca_dimensions
    ) +
      ggplot2::labs(
        title = paste0(args$project_name, ' - Elbow Plot'),
        subtitle = paste0('Showing first ', args$pca_dimensions, ' PCs')
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      )
    
    # Combine plots
    combined_plot <- (p1 | p2) +
      patchwork::plot_annotation(
        title = paste0(args$project_name, ' - PCA Analysis'),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16))
      )
    
    # Save plot
    plot_file <- file.path(output_dir, paste0(args$project_name, '_PCA_Summary.png'))
    ggplot2::ggsave(
      filename = plot_file,
      plot = combined_plot,
      width = 14,
      height = 7,
      units = 'in',
      dpi = config$project$dpi
    )
    
    logger_info(paste("PCA plots saved to:", plot_file))
    
  }, error = function(e) {
    logger_warning(paste("Error generating PCA plots:", e$message))
  })
}

# Function to generate UMAP plot
generate_umap_plot <- function(seurat_object, args, config, output_dir) {
  logger_info("Generating UMAP plot...")
  
  tryCatch({
    p <- Seurat::DimPlot(
      object = seurat_object,
      reduction = 'umap',
      label = TRUE,
      label.color = 'black',
      repel = TRUE,
      pt.size = config$visualization$point_size
    ) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = paste0(args$project_name, ' - UMAP Visualization'),
        subtitle = paste0('Clustering: ', args$clustering_algorithm, ' (Res: ', args$clustering_resolution, ')'),
        x = 'UMAP 1',
        y = 'UMAP 2'
      ) +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      )
    
    # Save plot
    plot_file <- file.path(output_dir, paste0(args$project_name, '_UMAP.png'))
    ggplot2::ggsave(
      filename = plot_file,
      plot = p,
      width = 10,
      height = 8,
      units = 'in',
      dpi = config$project$dpi
    )
    
    logger_info(paste("UMAP plot saved to:", plot_file))
    
  }, error = function(e) {
    logger_warning(paste("Error generating UMAP plot:", e$message))
  })
}

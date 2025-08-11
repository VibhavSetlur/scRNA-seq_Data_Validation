# Visualization Module for snRNA-seq Pipeline
# Handles additional visualizations and analysis

# Function to generate visualizations
generate_visualizations <- function(seurat_object, args, config, output_dir) {
  logger_info("Generating additional visualizations...")
  
  # Generate cluster statistics
  generate_cluster_statistics(seurat_object, args, config, output_dir)
  
  # Generate feature plots for top markers
  if (args$find_markers) {
    generate_feature_plots(seurat_object, args, config, output_dir)
  }
  
  # Generate heatmap
  generate_heatmap(seurat_object, args, config, output_dir)
  
  logger_success("Additional visualizations completed")
}

# Function to generate cluster statistics
generate_cluster_statistics <- function(seurat_object, args, config, output_dir) {
  logger_info("Generating cluster statistics...")
  
  tryCatch({
    # Get cluster information
    cluster_info <- table(Seurat::Idents(seurat_object))
    
    # Create cluster summary
    cluster_summary <- data.frame(
      Cluster = names(cluster_info),
      Cell_Count = as.numeric(cluster_info),
      Percent = round(as.numeric(cluster_info) / sum(cluster_info) * 100, 2)
    )
    
    # Save cluster summary
    summary_file <- file.path(output_dir, paste0(args$project_name, '_cluster_summary.csv'))
    write.csv(cluster_summary, file = summary_file, row.names = FALSE)
    
    # Create bar plot
    p <- ggplot2::ggplot(cluster_summary, ggplot2::aes(x = Cluster, y = Cell_Count, fill = Cluster)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes(label = Cell_Count), vjust = -0.5, size = 3) +
      ggplot2::labs(
        title = paste0(args$project_name, ' - Cell Counts per Cluster'),
        x = "Cluster",
        y = "Number of Cells"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    # Save plot
    plot_file <- file.path(output_dir, paste0(args$project_name, '_cluster_counts.png'))
    ggplot2::ggsave(
      filename = plot_file,
      plot = p,
      width = 10,
      height = 6,
      units = 'in',
      dpi = config$project$dpi
    )
    
    logger_info(paste("Cluster statistics saved to:", summary_file))
    
  }, error = function(e) {
    logger_warning(paste("Error generating cluster statistics:", e$message))
  })
}

# Function to generate feature plots
generate_feature_plots <- function(seurat_object, args, config, output_dir) {
  logger_info("Generating feature plots for top markers...")
  
  tryCatch({
    # Load markers if available
    markers_file <- file.path(output_dir, paste0(args$project_name, '_markers.csv'))
    
    if (!file.exists(markers_file)) {
      logger_warning("Markers file not found, skipping feature plots")
      return()
    }
    
    markers <- read.csv(markers_file)
    
    # Get top markers per cluster
    top_markers <- markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = 5, wt = avg_log2FC) %>%
      dplyr::pull(gene) %>%
      unique()
    
    # Limit to first 20 markers for visualization
    top_markers <- head(top_markers, 20)
    
    if (length(top_markers) > 0) {
      # Create feature plot
      p <- Seurat::FeaturePlot(
        object = seurat_object,
        features = top_markers,
        ncol = 4,
        pt.size = config$visualization$point_size
      ) +
        ggplot2::theme_classic() +
        ggplot2::labs(title = paste0(args$project_name, ' - Top Marker Genes'))
      
      # Save plot
      plot_file <- file.path(output_dir, paste0(args$project_name, '_feature_plots.png'))
      ggplot2::ggsave(
        filename = plot_file,
        plot = p,
        width = 16,
        height = 12,
        units = 'in',
        dpi = config$project$dpi
      )
      
      logger_info(paste("Feature plots saved to:", plot_file))
    }
    
  }, error = function(e) {
    logger_warning(paste("Error generating feature plots:", e$message))
  })
}

# Function to generate heatmap
generate_heatmap <- function(seurat_object, args, config, output_dir) {
  logger_info("Generating heatmap...")
  
  tryCatch({
    # Check if markers are available
    markers_file <- file.path(output_dir, paste0(args$project_name, '_markers.csv'))
    
    if (!file.exists(markers_file)) {
      logger_warning("Markers file not found, skipping heatmap")
      return()
    }
    
    markers <- read.csv(markers_file)
    
    # Get top markers per cluster
    top_markers <- markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = 10, wt = avg_log2FC) %>%
      dplyr::pull(gene) %>%
      unique()
    
    # Limit to first 50 markers for heatmap
    top_markers <- head(top_markers, 50)
    
    if (length(top_markers) > 0) {
      # Create heatmap
      p <- Seurat::DoHeatmap(
        object = seurat_object,
        features = top_markers,
        size = 3,
        angle = 45
      ) +
        ggplot2::labs(title = paste0(args$project_name, ' - Top Marker Genes Heatmap')) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14))
      
      # Save plot
      plot_file <- file.path(output_dir, paste0(args$project_name, '_heatmap.png'))
      ggplot2::ggsave(
        filename = plot_file,
        plot = p,
        width = 12,
        height = 10,
        units = 'in',
        dpi = config$project$dpi
      )
      
      logger_info(paste("Heatmap saved to:", plot_file))
    }
    
  }, error = function(e) {
    logger_warning(paste("Error generating heatmap:", e$message))
  })
}

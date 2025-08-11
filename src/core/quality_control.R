# Quality Control Module for snRNA-seq Pipeline
# Handles quality control filtering and visualization

# Function to perform quality control
perform_quality_control <- function(seurat_object, args, config, output_dir) {
  logger_info("Starting quality control...")
  
  # Calculate mitochondrial percentage
  seurat_object <- calculate_mt_percentage(seurat_object)
  
  # Get pre-filter statistics
  pre_filter_stats <- get_qc_statistics(seurat_object)
  
  # Generate pre-filter plots
  generate_qc_plots(seurat_object, args, config, output_dir, "pre_filter")
  
  # Filter doublets if SouporCell data provided
  if (!is.null(args$soupor_cell_doublet_input)) {
    seurat_object <- filter_doublets(seurat_object, args$soupor_cell_doublet_input)
  }
  
  # Apply QC filters
  seurat_object <- apply_qc_filters(seurat_object, args, config)
  
  # Get post-filter statistics
  post_filter_stats <- get_qc_statistics(seurat_object)
  
  # Generate post-filter plots
  generate_qc_plots(seurat_object, args, config, output_dir, "post_filter")
  
  # Generate QC summary
  generate_qc_summary(pre_filter_stats, post_filter_stats, output_dir, args$project_name)
  
  logger_success("Quality control completed")
  
  return(seurat_object)
}

# Function to calculate mitochondrial percentage
calculate_mt_percentage <- function(seurat_object) {
  logger_info("Calculating mitochondrial percentage...")
  
  tryCatch({
    seurat_object[['percent.mt']] <- Seurat::PercentageFeatureSet(
      object = seurat_object,
      pattern = '^mt-'
    )
    
    logger_success("Mitochondrial percentage calculated")
    return(seurat_object)
    
  }, error = function(e) {
    logger_warning(paste("Error calculating mitochondrial percentage:", e$message))
    # Add zero values if calculation fails
    seurat_object[['percent.mt']] <- rep(0, ncol(seurat_object))
    return(seurat_object)
  })
}

# Function to get QC statistics
get_qc_statistics <- function(seurat_object) {
  stats <- list(
    n_cells = ncol(seurat_object),
    n_genes = nrow(seurat_object),
    total_counts = sum(seurat_object$nCount_RNA),
    n_features_stats = list(
      min = min(seurat_object$nFeature_RNA, na.rm = TRUE),
      median = median(seurat_object$nFeature_RNA, na.rm = TRUE),
      max = max(seurat_object$nFeature_RNA, na.rm = TRUE)
    ),
    n_counts_stats = list(
      min = min(seurat_object$nCount_RNA, na.rm = TRUE),
      median = median(seurat_object$nCount_RNA, na.rm = TRUE),
      max = max(seurat_object$nCount_RNA, na.rm = TRUE)
    ),
    mt_percent_stats = list(
      min = min(seurat_object$percent.mt, na.rm = TRUE),
      median = median(seurat_object$percent.mt, na.rm = TRUE),
      max = max(seurat_object$percent.mt, na.rm = TRUE)
    )
  )
  
  return(stats)
}

# Function to filter doublets using SouporCell
filter_doublets <- function(seurat_object, souporcell_file) {
  logger_info("Filtering doublets using SouporCell data...")
  
  tryCatch({
    # Load SouporCell data
    souporcell_data <- load_souporcell_data(souporcell_file)
    
    if (is.null(souporcell_data)) {
      logger_warning("SouporCell data not available, skipping doublet filtering")
      return(seurat_object)
    }
    
    # Get barcodes of non-doublet cells
    non_doublets <- souporcell_data %>%
      dplyr::filter(status == 'singlet')
    
    # Find intersection with Seurat object
    valid_non_doublets <- intersect(non_doublets$barcode, colnames(seurat_object))
    
    if (length(valid_non_doublets) > 0) {
      # Record number of cells before filtering
      n_cells_before <- ncol(seurat_object)
      
      # Filter cells
      seurat_object <- subset(seurat_object, cells = valid_non_doublets)
      
      n_cells_after <- ncol(seurat_object)
      n_doublets_removed <- n_cells_before - n_cells_after
      
      logger_success(paste("Removed", n_doublets_removed, "doublets. Remaining cells:", n_cells_after))
    } else {
      logger_warning("No valid non-doublet cells found in SouporCell data")
    }
    
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error filtering doublets:", e$message))
    return(seurat_object)
  })
}

# Function to apply QC filters
apply_qc_filters <- function(seurat_object, args, config) {
  logger_info("Applying quality control filters...")
  
  # Get filter parameters
  min_features <- args$min_features
  min_counts <- args$min_counts
  max_features <- config$qc$max_features
  max_counts <- config$qc$max_counts
  max_mt_percent <- config$qc$max_mt_percent
  
  # Record initial cell count
  n_cells_before <- ncol(seurat_object)
  
  # Apply filters
  seurat_object <- subset(
    seurat_object,
    subset = nFeature_RNA >= min_features &
             nFeature_RNA <= max_features &
             nCount_RNA >= min_counts &
             nCount_RNA <= max_counts &
             percent.mt <= max_mt_percent
  )
  
  # Record final cell count
  n_cells_after <- ncol(seurat_object)
  n_cells_removed <- n_cells_before - n_cells_after
  
  logger_info(paste("QC filtering completed:"))
  logger_info(paste("  Cells before filtering:", n_cells_before))
  logger_info(paste("  Cells after filtering:", n_cells_after))
  logger_info(paste("  Cells removed:", n_cells_removed))
  logger_info(paste("  Filter parameters:"))
  logger_info(paste("    Min features:", min_features))
  logger_info(paste("    Max features:", max_features))
  logger_info(paste("    Min counts:", min_counts))
  logger_info(paste("    Max counts:", max_counts))
  logger_info(paste("    Max MT %:", max_mt_percent))
  
  return(seurat_object)
}

# Function to generate QC plots
generate_qc_plots <- function(seurat_object, args, config, output_dir, filter_type) {
  logger_info(paste("Generating", filter_type, "QC plots..."))
  
  project_name <- args$project_name
  
  # Create violin plots
  generate_qc_violin_plots(seurat_object, args, config, output_dir, filter_type)
  
  # Create feature scatter plots
  generate_qc_scatter_plots(seurat_object, args, config, output_dir, filter_type)
  
  # Create histogram plots
  generate_qc_histogram_plots(seurat_object, args, config, output_dir, filter_type)
  
  logger_success(paste(filter_type, "QC plots generated"))
}

# Function to generate violin plots
generate_qc_violin_plots <- function(seurat_object, args, config, output_dir, filter_type) {
  # Get statistics
  stats <- get_qc_statistics(seurat_object)
  
  # Create violin plots
  p1 <- Seurat::VlnPlot(
    object = seurat_object,
    features = 'nFeature_RNA',
    pt.size = 0
  ) +
    ggplot2::labs(
      title = 'Features per Cell',
      subtitle = paste0("Min: ", stats$n_features_stats$min, 
                       ", Median: ", round(stats$n_features_stats$median, 1),
                       ", Max: ", stats$n_features_stats$max)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8)
    ) +
    ggplot2::stat_summary(fun = "median", geom = "point", color = "black", size = 2)
  
  p2 <- Seurat::VlnPlot(
    object = seurat_object,
    features = 'nCount_RNA',
    pt.size = 0
  ) +
    ggplot2::labs(
      title = 'Counts per Cell',
      subtitle = paste0("Min: ", stats$n_counts_stats$min,
                       ", Median: ", round(stats$n_counts_stats$median, 1),
                       ", Max: ", stats$n_counts_stats$max)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8)
    ) +
    ggplot2::stat_summary(fun = "median", geom = "point", color = "black", size = 2)
  
  p3 <- Seurat::VlnPlot(
    object = seurat_object,
    features = 'percent.mt',
    pt.size = 0
  ) +
    ggplot2::labs(
      title = 'Percent Mitochondrial Reads',
      subtitle = paste0("Min: ", round(stats$mt_percent_stats$min, 1),
                       ", Median: ", round(stats$mt_percent_stats$median, 1),
                       ", Max: ", round(stats$mt_percent_stats$max, 1))
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8)
    ) +
    ggplot2::stat_summary(fun = "median", geom = "point", color = "black", size = 2)
  
  # Combine plots
  combined_plot <- (p1 | p2 | p3) +
    patchwork::plot_annotation(
      title = paste0(args$project_name, ' - ', stringr::str_to_title(filter_type), ' QC Metrics Distributions'),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16))
    )
  
  # Save plot
  plot_file <- file.path(output_dir, paste0(args$project_name, '_QC_VlnPlots_', stringr::str_to_title(filter_type), '_Combined.png'))
  ggplot2::ggsave(
    filename = plot_file,
    plot = combined_plot,
    width = 15,
    height = 6,
    units = 'in',
    dpi = config$project$dpi
  )
  
  logger_info(paste("Violin plots saved to:", plot_file))
}

# Function to generate scatter plots
generate_qc_scatter_plots <- function(seurat_object, args, config, output_dir, filter_type) {
  # Create feature scatter plots
  p1 <- Seurat::FeatureScatter(
    object = seurat_object,
    feature1 = 'nCount_RNA',
    feature2 = 'nFeature_RNA',
    pt.size = 0.5
  ) +
    ggplot2::labs(title = 'Features vs. Counts') +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12)) +
    ggplot2::scale_x_continuous(trans = 'log10')
  
  p2 <- Seurat::FeatureScatter(
    object = seurat_object,
    feature1 = 'nCount_RNA',
    feature2 = 'percent.mt',
    pt.size = 0.5
  ) +
    ggplot2::labs(title = 'Mitochondrial Percentage vs. Counts') +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12)) +
    ggplot2::scale_x_continuous(trans = 'log10')
  
  # Combine plots
  combined_plot <- (p1 | p2) +
    patchwork::plot_annotation(
      title = paste0(args$project_name, ' - ', stringr::str_to_title(filter_type), ' QC Feature Scatter (log10 Counts)'),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16))
    )
  
  # Save plot
  plot_file <- file.path(output_dir, paste0(args$project_name, '_QC_FeatureScatter_', stringr::str_to_title(filter_type), '_Combined.png'))
  ggplot2::ggsave(
    filename = plot_file,
    plot = combined_plot,
    width = 12,
    height = 6,
    units = 'in',
    dpi = config$project$dpi
  )
  
  logger_info(paste("Scatter plots saved to:", plot_file))
}

# Function to generate histogram plots
generate_qc_histogram_plots <- function(seurat_object, args, config, output_dir, filter_type) {
  # Create histogram data
  umi_data <- data.frame(UMI = seurat_object$nCount_RNA)
  stats <- get_qc_statistics(seurat_object)
  
  # Create histogram
  p <- ggplot2::ggplot(umi_data, ggplot2::aes(x = UMI)) +
    ggplot2::geom_histogram(
      fill = ifelse(filter_type == "pre_filter", '#4A76A8', '#356135'),
      color = 'dodgerblue',
      binwidth = 500
    ) +
    ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    ggplot2::labs(
      title = paste0(args$project_name, ' - UMI Counts per Cell Distribution (', stringr::str_to_title(filter_type), ')'),
      subtitle = paste0("Min: ", stats$n_counts_stats$min,
                       ", Median: ", round(stats$n_counts_stats$median, 1),
                       ", Max: ", stats$n_counts_stats$max),
      x = 'Total UMI Counts per Cell',
      y = 'Number of Cells'
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
    )
  
  # Save plot
  plot_file <- file.path(output_dir, paste0(args$project_name, '_QC_UMI_Histogram_', stringr::str_to_title(filter_type), '.png'))
  ggplot2::ggsave(
    filename = plot_file,
    plot = p,
    width = 10,
    height = 8,
    units = 'in',
    dpi = config$project$dpi
  )
  
  logger_info(paste("Histogram plot saved to:", plot_file))
}

# Function to generate QC summary
generate_qc_summary <- function(pre_filter_stats, post_filter_stats, output_dir, project_name) {
  logger_info("Generating QC summary...")
  
  # Create summary data
  summary_data <- data.frame(
    Metric = c("Number of Cells", "Number of Genes", "Total Counts"),
    Before_Filtering = c(
      pre_filter_stats$n_cells,
      pre_filter_stats$n_genes,
      pre_filter_stats$total_counts
    ),
    After_Filtering = c(
      post_filter_stats$n_cells,
      post_filter_stats$n_genes,
      post_filter_stats$total_counts
    )
  )
  
  # Calculate differences
  summary_data$Difference <- summary_data$Before_Filtering - summary_data$After_Filtering
  summary_data$Percent_Removed <- round((summary_data$Difference / summary_data$Before_Filtering) * 100, 2)
  
  # Save summary
  summary_file <- file.path(output_dir, paste0(project_name, '_QC_Summary.csv'))
  write.csv(summary_data, file = summary_file, row.names = FALSE)
  
  # Create bar plots
  generate_qc_bar_plots(summary_data, output_dir, project_name)
  
  logger_info(paste("QC summary saved to:", summary_file))
}

# Function to generate QC bar plots
generate_qc_bar_plots <- function(summary_data, output_dir, project_name) {
  # Cells bar plot
  cells_data <- data.frame(
    Filtering = factor(c("Before Filtering", "After Filtering"),
                      levels = c("Before Filtering", "After Filtering")),
    Count = c(summary_data$Before_Filtering[1], summary_data$After_Filtering[1])
  )
  
  p1 <- ggplot2::ggplot(cells_data, ggplot2::aes(x = Filtering, y = Count, fill = Filtering)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.5, color = 'black', size = 4) +
    ggplot2::labs(
      title = paste0(project_name, ' - Number of Cells Before and After Filtering'),
      x = "Filtering Status",
      y = "Number of Cells"
    ) +
    ggplot2::scale_fill_manual(values = c("Before Filtering" = "lightblue", "After Filtering" = "steelblue")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, size = 14))
  
  # Counts bar plot
  counts_data <- data.frame(
    Filtering = factor(c("Before Filtering", "After Filtering"),
                      levels = c("Before Filtering", "After Filtering")),
    Count = c(summary_data$Before_Filtering[3], summary_data$After_Filtering[3])
  )
  
  p2 <- ggplot2::ggplot(counts_data, ggplot2::aes(x = Filtering, y = Count, fill = Filtering)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.5, color = 'black', size = 4) +
    ggplot2::labs(
      title = paste0(project_name, ' - Total Reads Before and After Filtering'),
      x = "Filtering Status",
      y = "Total Reads"
    ) +
    ggplot2::scale_fill_manual(values = c("Before Filtering" = "lightgreen", "After Filtering" = "forestgreen")) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.1)),
      labels = scales::unit_format(unit = "M", scale = 1e-6)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, size = 14))
  
  # Save plots
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(project_name, '_QC_cells_summary.png')),
    plot = p1,
    width = 7,
    height = 6,
    units = 'in',
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(project_name, '_QC_reads_summary.png')),
    plot = p2,
    width = 7,
    height = 6,
    units = 'in',
    dpi = 300
  )
  
  logger_info("QC bar plots generated")
}

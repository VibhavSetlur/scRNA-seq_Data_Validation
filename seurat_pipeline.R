#!/usr/bin/env Rscript

# Seurat Single-Cell RNA-seq Pipeline
# Description: This script performs a standard snRNA-cell analysis pipeline using Seurat.
#              It can take either a Seurat RDS file or a raw H5 data file as input.
#              If an H5 file is provided, QC filtering thresholds can be specified
#              via command-line arguments.

# Check if necessary libraries are installed and install them if not
list_of_packages = c("Seurat", "tidyverse", "patchwork", "argparse", "plotly", "htmlwidgets", "scales", "ggrepel")

cat("Checking and installing necessary packages...\n")
for (package_name in list_of_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    cat(paste("Installing package:", package_name, "...\n"))
    install.packages(package_name)
    if (!requireNamespace(package_name, quietly = TRUE)) {
      stop(paste("Error: Package", package_name, "could not be installed."))
    }
  } else {
    cat(paste("Package", package_name, "is already installed.\n"))
  }
}

cat("All necessary packages should now be installed and loaded.\n")

# Load necessary libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(argparse)
library(plotly)
library(htmlwidgets)
library(scales)
library(ggrepel)


# Function to calculate summary statistics
calculate_summary_stats <- function(data_vector) {
  # Check if data_vector is empty or contains only NA values
  if (length(data_vector) == 0 || all(is.na(data_vector))) {
    return(list(min = NA, median = NA, max = NA))
  }
  summary_stats <- list(
    min = min(data_vector, na.rm = TRUE),
    median = median(data_vector, na.rm = TRUE),
    max = max(data_vector, na.rm = TRUE)
  )
  return(summary_stats)
}


# Set up argument parsing
parser = ArgumentParser(description = 'Run Seurat snRNA-cell pipeline.')

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
parser$add_argument('--soupor_cell_doublet_input',
                    type = 'character',
                    help = 'SouporCell clusters output to identify and remove doublets from the data.')
parser$add_argument('--min_features',
                    type = 'integer',
                    default = 200,
                    help = 'Minimum number of features (genes) per cell for QC filtering (only if H5 input). Default value = 200')
parser$add_argument('--n_variable_features',
                    type = 'integer',
                    default = 2000,
                    help = 'Number of variable features to find. Default value = 2000')
parser$add_argument('--normalization_method',
                    type = 'character',
                    default = 'CLR',
                    choices = c('LogNormalize', 'RC', 'CLR'),
                    help = 'Type of normalization method to be performed (LogNormalize, RC, CLR). Default value = CLR')
parser$add_argument('--min_counts',
                    type = 'integer',
                    default = 1000,
                    help = 'Minimum number of counts per cell for QC filtering (only if H5 input). Default value = 1000')
parser$add_argument('--scaling_method',
                    type = 'character',
                    default = 'negbinom',
                    choices = c('negbinom', 'linear'),
                    help = 'Model to use for scaling data (negbinom or linear). negbinom is recommended for UMI datasets. Default value = negbinom')
parser$add_argument('--pca_dimensions',
                    type = 'integer',
                    default = 15,
                    help = 'Number of principal components to compute. Default value = 15')
parser$add_argument('--clustering_resolution',
                    type = 'double',
                    default = 0.5,
                    help = 'Resolution parameter for FindClusters, influences the number of clusters. Default value = 0.5')
parser$add_argument('--clustering_algorithm',
                    type = 'character',
                    default = 'leiden',
                    choices = c('louvain', 'multilevel', 'leiden', 'slm'),
                    help = 'Clustering algorithm to use (louvain, multilevel, leiden, slm). Leiden generally optimizes modularity better for larger datasets. Default value = leiden')
parser$add_argument('--working_dir',
                    type = 'character',
                    default = '.',
                    help = 'Working directory where output files will be saved.')
parser$add_argument('--find_markers',
                    type = 'logical',
                    default = TRUE,
                    help = 'Boolean value to determine whether to find cluster markers. Default value = TRUE')

args = parser$parse_args()

# Set working directory
setwd(args$working_dir)

# Load data
seurat_object = NULL
input_type = NULL

cat("Loading data...\n")
if (!is.null(args$rds_input) && file.exists(args$rds_input)) {
  
  cat('Loading Seurat object from:', args$rds_input, '\n')
  seurat_object = readRDS(file = args$rds_input)
  input_type = 'rds'
  cat("Seurat object loaded successfully.\n")
  
} else if (!is.null(args$h5_input) && file.exists(args$h5_input)) {
  
  cat('Loading raw data from:', args$h5_input, '\n')
  h5_data = Read10X_h5(filename = args$h5_input)
  seurat_object = CreateSeuratObject(counts = h5_data,
                                     project = args$project_name)
  input_type = 'h5'
  cat("Raw data loaded and Seurat object created.\n")
  
} else {
  
  stop("Please provide a valid existing RDS input file or H5 input file.")
  
}

#---------------------------------------------------------------------------------------
# Quality control (performed if Seurat object is loaded or created)
if (!is.null(seurat_object)) {
  cat('Performing initial quality control...\n')
  
  cat("Calculating mitochondrial percentage...\n")
  seurat_object[['percent.mt']] = PercentageFeatureSet(object = seurat_object,
                                                       pattern = '^mt-')
  cat("Mitochondrial percentage calculated.\n")
  
  # Get the number of cells and reads before filtering
  n_cells_before = ncol(seurat_object)
  n_reads_before = sum(seurat_object$nCount_RNA)
  
  # Calculate pre-filter summary statistics
  n_features_pre_stats <- calculate_summary_stats(seurat_object$nFeature_RNA)
  n_counts_pre_stats <- calculate_summary_stats(seurat_object$nCount_RNA)
  percent_mt_pre_stats <- calculate_summary_stats(seurat_object$percent.mt)
  
  # Format pre-filter summary statistics for subtitle
  n_features_pre_text <- paste0("Min: ", n_features_pre_stats$min, ", Median: ", round(n_features_pre_stats$median, 1), ", Max: ", n_features_pre_stats$max)
  n_counts_pre_text <- paste0("Min: ", n_counts_pre_stats$min, ", Median: ", round(n_counts_pre_stats$median, 1), ", Max: ", n_counts_pre_stats$max)
  percent_mt_pre_text <- paste0("Min: ", round(percent_mt_pre_stats$min, 1), ", Median: ", round(percent_mt_pre_stats$median, 1), ", Max: ", round(percent_mt_pre_stats$max, 1))
  
  
  cat('Generating QC VlnPlots (Static, Pre-filter)...\n')
  # View distribution of feature and counts per cell - Separate VlnPlots (Static)
  p1_nFeature_pre = VlnPlot(object = seurat_object, features = 'nFeature_RNA', pt.size = 0) + # Increased pt.size slightly for visibility if not too many points
    labs(title = 'Features per Cell', subtitle = n_features_pre_text) + # Added stats to subtitle
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12), # Center and size title for combined view
          plot.subtitle = element_text(hjust = 0.5, size = 8)) + # Center and size subtitle
    stat_summary(fun = "median", geom = "point", color = "black", size = 2) # Add median point
  
  p1_nCount_pre = VlnPlot(object = seurat_object, features = 'nCount_RNA', pt.size = 0) +
    labs(title = 'Counts per Cell', subtitle = n_counts_pre_text) + # Added stats to subtitle
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12), # Center and size title for combined view
          plot.subtitle = element_text(hjust = 0.5, size = 8)) + # Center and size subtitle
    stat_summary(fun = "median", geom = "point", color = "black", size = 2) # Add median point
  
  p1_percentmt_pre = VlnPlot(object = seurat_object, features = 'percent.mt', pt.size = 0) +
    labs(title = 'Percent Mitochondrial Reads', subtitle = percent_mt_pre_text) + # Added stats to subtitle
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12), # Center and size title for combined view
          plot.subtitle = element_text(hjust = 0.5, size = 8)) + # Center and size subtitle
    stat_summary(fun = "median", geom = "point", color = "black", size = 2) # Add median point
  
  # Combine Pre-filter VlnPlots using patchwork
  combined_vln_pre = (p1_nFeature_pre | p1_nCount_pre | p1_percentmt_pre) +
    plot_annotation(title = paste0(args$project_name, ' - Pre-filter QC Metrics Distributions'),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) # Add overall title
  
  
  # Save combined Pre-filter VlnPlots
  ggsave(filename = paste0(args$project_name, '_QC_VlnPlots_PreFilter_Combined.png'), plot = combined_vln_pre, width = 15, height = 6, units = 'in', dpi = 300)
  cat('Combined QC VlnPlots (static, pre-filter) saved.\n')
  
  
  cat('Generating QC FeatureScatter plots (Static, Pre-filter)...\n')
  # View number of counts by number of features - FeatureScatter (Static)
  p2_pre = FeatureScatter(object = seurat_object,
                          feature1 = 'nCount_RNA',
                          feature2 = 'nFeature_RNA',
                          pt.size = 0.5) +
    labs(title = 'Features vs. Counts') + # Simplified title for combined plot
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) + # Center and size title for combined view
    scale_x_continuous(trans = 'log10') # Apply log10 scale to x-axis
  
  p2_mt_pre = FeatureScatter(object = seurat_object,
                             feature1 = 'nCount_RNA',
                             feature2 = 'percent.mt',
                             pt.size = 0.5) +
    labs(title = 'Mitochondrial Percentage vs. Counts') + # Simplified title for combined plot
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    scale_x_continuous(trans = 'log10') # Apply log10 scale to x-axis
  
  
  # Combine Pre-filter FeatureScatter plots using patchwork
  combined_scatter_pre = (p2_pre | p2_mt_pre) +
    plot_annotation(title = paste0(args$project_name, ' - Pre-filter QC Feature Scatter (log10 Counts)'), # Indicate log scale
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) # Add overall title
  
  # Save combined Pre-filter FeatureScatter plots
  ggsave(filename = paste0(args$project_name, '_QC_FeatureScatter_PreFilter_Combined.png'), plot = combined_scatter_pre, width = 12, height = 6, units = 'in', dpi = 300)
  cat('Combined QC FeatureScatter plots (static, pre-filter, log10 counts) saved.\n')
  
  
  cat('Generating QC Histogram (Pre-filter)...\n')
  # View cell distribution (UMI counts per cell) - Histogram (Pre-filter) (Static)
  p3_data <- data.frame(UMI = seurat_object$nCount_RNA)
  umi_pre_stats <- calculate_summary_stats(p3_data$UMI)
  umi_pre_text <- paste0("Min: ", umi_pre_stats$min, ", Median: ", round(umi_pre_stats$median, 1), ", Max: ", umi_pre_stats$max)
  
  p3 = ggplot(p3_data, aes(x = UMI)) +
    geom_histogram(fill = '#4A76A8', color = 'black', binwidth = 500) + # Changed fill color
    scale_x_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) + # Keep linear scale with unit format
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    labs(title = paste0(args$project_name, ' - UMI Counts per Cell Distribution (Pre-filter)'),
         subtitle = umi_pre_text, # Added stats to subtitle
         x = 'Total UMI Counts per Cell',
         y = 'Number of Cells') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  # Save Pre-filter Histogram
  ggsave(filename = paste0(args$project_name, '_QC_UMI_Histogram_PreFilter.png'), plot = p3, width = 7, height = 6, units = 'in', dpi = 300)
  cat('QC Histogram (Pre-filter, static) saved.\n')
  
  # Filter for doublets if souporcell input is provided
  if (!is.null(args$soupor_cell_doublet_input) && file.exists(args$soupor_cell_doublet_input)) {
    cat('Filtering doublets based on SouporCell output...\n')
    clusters = read.csv(args$soupor_cell_doublet_input,
                        sep = '\t')
    
    # Get barcodes of non-doublet cells
    if ("barcode" %in% colnames(clusters)) {
      non_doublets <- clusters %>%
        filter(status != 'doublet')
      
      valid_non_doublets <- intersect(non_doublets$barcode, colnames(seurat_object))
      
      if(length(valid_non_doublets) > 0) {
        # Record number of cells before doublet removal
        n_cells_before_doublet_filter = ncol(seurat_object)
        seurat_object <- subset(seurat_object, cells = valid_non_doublets)
        n_cells_after_doublet_filter = ncol(seurat_object)
        n_doublets_removed = n_cells_before_doublet_filter - n_cells_after_doublet_filter
        cat(paste("Removed", n_doublets_removed, "doublets based on SouporCell. Number of cells remaining:", ncol(seurat_object), "\n"))
        
      } else {
        cat("Warning: No SouporCell barcodes matched Seurat object cell names or all matched cells were doublets/filtered. Skipping doublet filtering.\n")
      }
      
    } else {
      cat("Error: SouporCell input file does not contain a 'barcode' column. Skipping doublet filtering.\n")
    }
  } else if (!is.null(args$soupor_cell_doublet_input) && !file.exists(args$soupor_cell_doublet_input)) {
    cat(paste("Warning: SouporCell doublet input file not found at", args$soupor_cell_doublet_input, ". Skipping doublet filtering.\n"))
  }
  
  
  # Filter cells based on user-provided thresholds (only if H5 input)
  if (input_type == 'h5') {
    min_features = args$min_features
    min_counts = args$min_counts
    
    cat(paste('Filtering cells with at least', min_features, 'features and', min_counts, 'counts...\n'))
    
    # Check if there are cells to filter before attempting subset
    if (ncol(seurat_object) > 0) {
      seurat_object = subset(x = seurat_object,
                             subset = nFeature_RNA >= min_features & nCount_RNA >= min_counts)
      cat(paste("Filtered by min_features (", min_features, ") and min_counts (", min_counts, "). Number of cells remaining:", ncol(seurat_object), "\n"))
    } else {
      cat("Warning: No cells remaining after previous filtering steps. Skipping min_features and min_counts filtering.\n")
    }
    
    
    # Only generate post-filter plots if cells remain
    if (ncol(seurat_object) > 0) {
      # Calculate post-filter summary statistics
      n_features_post_stats <- calculate_summary_stats(seurat_object$nFeature_RNA)
      n_counts_post_stats <- calculate_summary_stats(seurat_object$nCount_RNA)
      percent_mt_post_stats <- calculate_summary_stats(seurat_object$percent.mt)
      
      # Format post-filter summary statistics for subtitle
      n_features_post_text <- paste0("Min: ", n_features_post_stats$min, ", Median: ", round(n_features_post_stats$median, 1), ", Max: ", n_features_post_stats$max)
      n_counts_post_text <- paste0("Min: ", n_counts_post_stats$min, ", Median: ", round(n_counts_post_stats$median, 1), ", Max: ", n_counts_post_stats$max)
      percent_mt_post_text <- paste0("Min: ", round(percent_mt_post_stats$min, 1), ", Median: ", round(percent_mt_post_stats$median, 1), ", Max: ", round(percent_mt_post_stats$max, 1))
      
      
      cat('Generating QC VlnPlots (Static, Post-filter)...\n')
      p1_nFeature_post = VlnPlot(object = seurat_object, features = 'nFeature_RNA', pt.size = 0) +
        labs(title = paste0('Features per Cell'),
             subtitle = paste0('Min Features: ', min_features, '\n', n_features_post_text)) + # Added filter and stats info
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              plot.subtitle = element_text(hjust = 0.5, size = 8)) +
        stat_summary(fun = "median", geom = "point", color = "black", size = 2) # Add median point
      
      p1_nCount_post = VlnPlot(object = seurat_object, features = 'nCount_RNA', pt.size = 0) +
        labs(title = paste0('Counts per Cell'),
             subtitle = paste0('Min Counts: ', min_counts, '\n', n_counts_post_text)) + # Added filter and stats info
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              plot.subtitle = element_text(hjust = 0.5, size = 8)) +
        stat_summary(fun = "median", geom = "point", color = "black", size = 2) # Add median point
      
      p1_percentmt_post = VlnPlot(object = seurat_object, features = 'percent.mt', pt.size = 0) +
        labs(title = paste0('Percent Mitochondrial Reads'),
             subtitle = percent_mt_post_text) + # Added stats info (MT filter threshold not applied)
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              plot.subtitle = element_text(hjust = 0.5, size = 8)) +
        stat_summary(fun = "median", geom = "point", color = "black", size = 2) # Add median point
      
      # Combine Post-filter VlnPlots using patchwork
      combined_vln_post = (p1_nFeature_post | p1_nCount_post | p1_percentmt_post) +
        plot_annotation(title = paste0(args$project_name, ' - Post-filter QC Metrics Distributions'),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) # Add overall title
      
      ggsave(filename = paste0(args$project_name, '_QC_VlnPlots_PostFilter_Combined.png'), plot = combined_vln_post, width = 15, height = 6, units = 'in', dpi = 300)
      cat('Combined QC VlnPlots (static, post-filter) saved.\n')
      
      cat('Generating QC FeatureScatter plots (Static, Post-filter)...\n')
      p2_post = FeatureScatter(object = seurat_object,
                               feature1 = 'nCount_RNA',
                               feature2 = 'nFeature_RNA',
                               pt.size = 0.5) +
        labs(title = 'Features vs. Counts',
             subtitle = paste0('Min Features: ', min_features, ', Min Counts: ', min_counts)) + # Added filter info in subtitle
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              plot.subtitle = element_text(hjust = 0.5, size = 9)) +
        scale_x_continuous(trans = 'log10') # Apply log10 scale to x-axis
      
      p2_mt_post = FeatureScatter(object = seurat_object,
                                  feature1 = 'nCount_RNA',
                                  feature2 = 'percent.mt',
                                  pt.size = 0.5) +
        labs(title = 'Mitochondrial Percentage vs. Counts',
             subtitle = paste0('Min Counts: ', min_counts)) + # MT filter not applied, but counts filter was
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              plot.subtitle = element_text(hjust = 0.5, size = 9)) +
        scale_x_continuous(trans = 'log10') # Apply log10 scale to x-axis
      
      # Combine Post-filter FeatureScatter plots using patchwork
      combined_scatter_post = (p2_post | p2_mt_post) +
        plot_annotation(title = paste0(args$project_name, ' - Post-filter QC Feature Scatter (log10 Counts)'), # Indicate log scale
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) # Add overall title
      
      
      ggsave(filename = paste0(args$project_name, '_QC_FeatureScatter_PostFilter_Combined.png'), plot = combined_scatter_post, width = 12, height = 6, units = 'in', dpi = 300)
      cat('Combined QC FeatureScatter plots (static, post-filter, log10 counts) saved.\n')
      
      
      # Only generate Post-filter histogram if filtering was applied (i.e., input was h5)
      cat('Generating QC Histogram (Post-filter)...\n')
      p4_data <- data.frame(UMI = seurat_object$nCount_RNA)
      umi_post_stats <- calculate_summary_stats(p4_data$UMI)
      umi_post_text <- paste0("Min: ", umi_post_stats$min, ", Median: ", round(umi_post_stats$median, 1), ", Max: ", umi_post_stats$max)
      
      p4 = ggplot(p4_data, aes(x = UMI)) +
        geom_histogram(fill = '#356135', color = 'black', binwidth = 500) + # Changed fill color
        scale_x_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) + # Keep linear scale with unit format
        scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
        labs(title = paste0(args$project_name, ' - UMI Counts per Cell Distribution (Post-filter)'),
             subtitle = paste0('Filtered (Min Features: ', args$min_features, ', Min Counts: ', args$min_counts, ')\n', umi_post_text), # Added filter and stats info
             x = 'Total UMI Counts per Cell',
             y = 'Number of Cells') +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              plot.subtitle = element_text(hjust = 0.5, size = 10)) # Center title and subtitle
      
      # Save Post-filter Histogram
      ggsave(filename = paste0(args$project_name, '_QC_UMI_Histogram_PostFilter.png'), plot = p4, width = 7, height = 6, units = 'in', dpi = 300)
      cat('QC Histogram (Post-filter, static) saved.\n')
    } else {
      cat("Warning: No cells remaining after filtering by min_features and min_counts. Skipping post-filter QC plot generation.\n")
    }
    
  }
  
  
  # Get the number of cells and reads after all filtering steps
  n_cells_after = ncol(seurat_object)
  n_reads_after = if(ncol(seurat_object) > 0) sum(seurat_object$nCount_RNA) else 0 # Handle case where no cells remain
  
  # Create dataframes for cell and read counts
  plot_data_cells <- data.frame(
    Filtering = factor(c("Before Filtering", "After Filtering"),
                       levels = c("Before Filtering", "After Filtering")),
    Count = c(n_cells_before, n_cells_after)
  )
  
  plot_data_reads <- data.frame(
    Filtering = factor(c("Before Filtering", "After Filtering"),
                       levels = c("Before Filtering", "After Filtering")),
    Count = c(n_reads_before, n_reads_after)
  )
  
  cat('Generating QC count summary bar plots (Static)...\n')
  # Create the bar plots using ggplot2 - Enhanced Bar Plots (Static)
  cell_number_plot <- ggplot(plot_data_cells, aes(x = Filtering, y = Count, fill = Filtering)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5, color = 'black', size = 4) + # Adjusted vjust to be above bar, black color, size 4
    labs(title = paste0(args$project_name, ' - Number of Cells Before and After Filtering'),
         x = "Filtering Status",
         y = "Number of Cells") +
    scale_fill_manual(values = c("Before Filtering" = "lightblue", "After Filtering" = "steelblue")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Ensure space for text above bar
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
  
  read_number_plot <- ggplot(plot_data_reads, aes(x = Filtering, y = Count, fill = Filtering)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5, color = 'black', size = 4) + # Adjusted vjust to be above bar, black color, size 4
    labs(title = paste0(args$project_name, ' - Total Reads Before and After Filtering'),
         x = "Filtering Status",
         y = "Total Reads") +
    scale_fill_manual(values = c("Before Filtering" = "lightgreen", "After Filtering" = "forestgreen")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = scales::unit_format(unit = "M", scale = 1e-6)) + # Scale y-axis for reads
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
  
  # Save bar plots
  ggsave(filename = paste0(args$project_name, '_QC_cells_summary.png'), plot = cell_number_plot, width = 7, height = 6, units = 'in', dpi = 300)
  ggsave(filename = paste0(args$project_name, '_QC_reads_summary.png'), plot = read_number_plot, width = 7, height = 6, units = 'in', dpi = 300)
  cat('QC count summary bar plots (static) saved.\n')
  
}


# Continue with downstream analysis only if Seurat object is still valid and has cells
if (!is.null(seurat_object) && ncol(seurat_object) > 0) {
  
  #---------------------------------------------------------------------------------------
  
  cat(paste('Normalizing data using', args$normalization_method , '...\n'))
  cat("Starting NormalizeData...\n")
  # Normalize data
  normalization_method = args$normalization_method
  scale_factor_val = NULL
  
  if (normalization_method == 'LogNormalize' || normalization_method == 'RC') {
    cat('Using default scale factor for LogNormalize or RC.\n')
  } else if (normalization_method == 'CLR') {
    scale_factor_val = 10000 # Default scale factor for CLR in Seurat
  }
  
  seurat_object = NormalizeData(object = seurat_object,
                                normalization.method = normalization_method,
                                scale.factor = scale_factor_val,
                                margin = 2) # cells
  cat("NormalizeData complete.\n")
  
  # Identify variable features
  cat('Finding variable features...\n')
  cat("Starting FindVariableFeatures...\n")
  
  nfeatures = args$n_variable_features
  
  # Check if there are enough features to find variable features
  if (nrow(seurat_object) < nfeatures) {
    cat(paste0("Warning: Number of features in object (", nrow(seurat_object), ") is less than requested number of variable features (", nfeatures, "). Setting n_variable_features to total features.\n"))
    nfeatures = nrow(seurat_object)
    if (nfeatures == 0) {
      cat("Error: No features found in object. Cannot find variable features. Skipping downstream analysis.\n")
      seurat_object = NULL # Set object to NULL to skip subsequent steps
    }
  }
  
  if (!is.null(seurat_object)) {
    seurat_object = FindVariableFeatures(object = seurat_object,
                                         selection.method = 'vst',
                                         nfeatures = nfeatures)
    cat("FindVariableFeatures complete.\n")
    
    top_variable_features = head(VariableFeatures(seurat_object), 10)
    
    cat('Generating Variable Feature Plot (Interactive)...\n')
    # Variable Feature Plot (Interactive with Hover)
    # Create the base ggplot object
    p5_static <- tryCatch({
      VariableFeaturePlot(object = seurat_object) +
        theme_classic() +
        labs(title = paste0(args$project_name, ' - Variable Features (Top ', length(VariableFeatures(seurat_object)), ')'), # Use actual number of variable features found
             x = 'Average Expression',
             y = 'Variance (vst)') +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
    }, error = function(e) {
      cat(paste0("Error generating static Variable Feature Plot: ", e$message, "\n"))
      return(NULL) # Return NULL if plotting fails
    })
    
    
    # Add static labels for the top 10 features using LabelPoints and ggrepel
    if (!is.null(p5_static) && length(top_variable_features) > 0) {
      tryCatch({
        # Add labels using LabelPoints
        p5_static <- LabelPoints(
          plot = p5_static,
          points = top_variable_features,
          repel = TRUE
        )
      }, error = function(e) {
        cat(paste0("Error adding labels to static Variable Feature Plot: ", e$message, "\n"))
      })
    } else {
      cat("Warning: No variable features were found. Cannot add labels to Variable Feature Plot.\n")
    }
    
    
    
    # Convert the static plot to interactive plotly object with gene name tooltip
    if (!is.null(p5_static)) {
      cat('Generating Interactive Variable Feature Plot (Plotly)...\n')
      p5_interactive <- tryCatch({
        ggplotly(p5_static, tooltip = c("name", "x", "y")) %>%
          layout(title = paste0(args$project_name, ' - Variable Features (Top ', length(VariableFeatures(seurat_object)), ')'))
      }, error = function(e) {
        cat(paste0("Error converting Variable Feature Plot to interactive Plotly: ", e$message, "\n"))
        return(NULL) # Return NULL if plotly conversion fails
      })
      
      # Save interactive Variable Features plot as HTML
      if (!is.null(p5_interactive)) {
        saveWidget(p5_interactive,
                   file = paste0(args$project_name, '_Variable_Features.html'),
                   selfcontained = TRUE)
        cat('Variable Feature Plot (interactive HTML) saved.\n')
      } else {
        cat('Skipping saving interactive Variable Feature Plot due to previous error.\n')
      }
    } else {
      cat('Skipping generation and saving of Variable Feature Plot due to previous error.\n')
    }
    
  }
  
  
  #---------------------------------------------------------------------------------------
  # Continue with downstream analysis only if Seurat object is still valid and has cells and variable features
  if (!is.null(seurat_object) && ncol(seurat_object) > 0 && length(VariableFeatures(seurat_object)) > 0) {
    cat('Scaling data using', args$scaling_method, 'model...\n')
    cat("Starting ScaleData...\n")
    # Scale data
    seurat_object = ScaleData(object = seurat_object,
                              features = VariableFeatures(seurat_object), # Scale only variable features
                              do.scale = TRUE,
                              do.center = TRUE,
                              model.use = args$scaling_method,
                              assay = 'RNA')
    cat("ScaleData complete.\n")
    
    
    #---------------------------------------------------------------------------------------
    cat('Performing dimensionality reduction (PCA)....\n')
    cat("Starting RunPCA...\n")
    # Perform dimensionality reduction (PCA)
    pca_dimensions = args$pca_dimensions
    # Check if the number of variable features is less than the number of requested PCs
    if (length(VariableFeatures(seurat_object)) < pca_dimensions) {
      cat(paste0("Warning: Number of variable features (", length(VariableFeatures(seurat_object)), ") is less than requested PCA dimensions (", pca_dimensions, "). Setting PCA dimensions to number of variable features.\n"))
      pca_dimensions = length(VariableFeatures(seurat_object))
      if (pca_dimensions == 0) {
        cat("Error: No variable features found. Cannot perform PCA. Skipping downstream analysis.\n")
        seurat_object = NULL # Set object to NULL to skip subsequent steps
      }
    }
    
    if (!is.null(seurat_object)) {
      # Check if there are enough cells for PCA (should be > 1 cell)
      if (ncol(seurat_object) > 1) {
        seurat_object = RunPCA(object = seurat_object,
                               npcs = pca_dimensions,
                               rev.pca = TRUE,
                               weight.by.var = TRUE,
                               verbose = FALSE)
        cat("RunPCA complete.\n")
        
        cat('Generating PCA plots (Static)...\n')
        # PCA DimPlot (Static)
        p6 <- tryCatch({
          DimPlot(object = seurat_object, reduction = 'pca') +
            labs(title = paste0(args$project_name, ' - PCA Plot'), subtitle = paste0('Computed Dimensions: ', pca_dimensions)) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5, size = 14), plot.subtitle = element_text(hjust = 0.5, size = 10))
        }, error = function(e) {
          cat(paste0("Error generating PCA DimPlot: ", e$message, "\n"))
          return(NULL)
        })
        
        # ElbowPlot (Static)
        p7 <- tryCatch({
          ElbowPlot(object = seurat_object, reduction = 'pca', ndims = pca_dimensions) +
            labs(title = paste0(args$project_name, ' - Elbow Plot'), subtitle = paste0('Showing first ', pca_dimensions, ' PCs')) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5, size = 14), plot.subtitle = element_text(hjust = 0.5, size = 10))
        }, error = function(e) {
          cat(paste0("Error generating PCA ElbowPlot: ", e$message, "\n"))
          return(NULL)
        })
        
        
        # DimHeatmap (Static) - Wrapped in tryCatch
        cat("Attempting to generate DimHeatmap...\n")
        p8 <- tryCatch({
          # Check if PCA calculation was successful and if number of dims is valid
          if (!is.null(seurat_object@reductions$pca) &&
              inherits(seurat_object@reductions$pca, "DimReduc") &&
              !is.null(seurat_object@reductions$pca@feature.loadings) &&
              ncol(seurat_object@reductions$pca@feature.loadings) >= min(pca_dimensions, 15) &&
              min(pca_dimensions, 15) > 0) {
            
            DimHeatmap(object = seurat_object,
                       dims = 1:min(pca_dimensions, 15),
                       nfeatures = 30,
                       slot = 'scale.data',
                       fast = TRUE,
                       combine = TRUE,
                       balanced = TRUE) +
              labs(title = paste0(args$project_name, ' - Top Features per PC (Scaled Data)')) +
              theme(plot.title = element_text(hjust = 0.5, size = 14))
          } else {
            cat("Warning: PCA reduction object is invalid or insufficient dimensions for DimHeatmap.\n")
            return(NULL) # Explicitly return NULL if checks fail within tryCatch
          }
        }, error = function(e) {
          cat(paste0("Error generating PCA DimHeatmap plot: ", e$message, "\n"))
          cat("Skipping DimHeatmap plot generation.\n")
          return(NULL) # Return NULL if an error occurs
        })
        
        # Save DimHeatmap (static) if it was successfully created
        if (!is.null(p8)) {
          ggsave(filename = paste0(args$project_name, '_PCA_DimHeatmap.png'),
                 plot = p8,
                 width = 12,
                 height = 10,
                 units = 'in',
                 dpi = 300)
          cat('PCA DimHeatmap plot (static) saved.\n')
        } else {
          cat('Skipping saving PCA DimHeatmap plot due to previous error or warning.\n')
        }
        
        
        # Combine PCA DimPlot and ElbowPlot for saving using patchwork
        if (!is.null(p6) && !is.null(p7)) {
          combined_plot_2 = (p6 | p7)
          # Save combined PCA plots (static)
          ggsave(filename = paste0(args$project_name, '_PCA_Summary.png'),
                 plot = combined_plot_2,
                 width = 14,
                 height = 7,
                 units = 'in',
                 dpi = 300)
          cat('Combined PCA plots (static) saved.\n')
        } else {
          cat('Skipping saving combined PCA plots due to previous errors.\n')
        }
        
      } else {
        cat("Error: Not enough cells remaining after filtering to perform PCA (need > 1 cell). Skipping downstream analysis.\n")
        seurat_object = NULL
      }
      
    }
  } else {
    cat("Skipping scaling and PCA steps as Seurat object is not valid, has no cells, or no variable features.\n")
  }
  
  
  #---------------------------------------------------------------------------------------
  # Continue with downstream analysis only if Seurat object is still valid and has cells and PCA reduction
  # Check if PCA reduction exists and has valid dimensions before proceeding with clustering and UMAP
  pca_is_valid <- !is.null(seurat_object) &&
    ncol(seurat_object) > 0 &&
    !is.null(seurat_object@reductions$pca) &&
    inherits(seurat_object@reductions$pca, "DimReduc") &&
    !is.null(seurat_object@reductions$pca@feature.loadings) &&
    ncol(seurat_object@reductions$pca@feature.loadings) > 0
  
  if (pca_is_valid) {
    cat('Performing clustering and UMAP...\n')
    cat("Starting FindNeighbors...\n")
    # Perform clustering and UMAP
    # Ensure pca_dims is valid based on the actual number of computed PCs
    valid_pca_dims = 1:min(args$pca_dimensions, ncol(seurat_object@reductions$pca@feature.loadings))
    
    # Check if there are valid PCA dimensions to use for FindNeighbors
    if (length(valid_pca_dims) > 0 && max(valid_pca_dims) > 0) {
      seurat_object = FindNeighbors(object = seurat_object,
                                    dims = valid_pca_dims,
                                    nn.method = 'annoy',
                                    verbose = FALSE)
      cat("FindNeighbors complete.\n")
      
      
      # Save Neighbors Graph (Optional, mainly for advanced users/debugging)
      cat('Saving Neighbors graphs...\n')
      if (!is.null(seurat_object@graphs$RNA_nn)) {
        saveRDS(seurat_object@graphs$RNA_nn, file = paste0(args$project_name, '_neighbors_graph.rds'))
        cat('Neighbors graph saved.\n')
      } else {
        cat("Warning: RNA_nn graph not found. Skipping saving neighbors graph.\n")
      }
      if (!is.null(seurat_object@graphs$RNA_snn)) {
        saveRDS(seurat_object@graphs$RNA_snn, file = paste0(args$project_name, '_shared_nearest_neighbors_graph.rds'))
        cat('Shared nearest neighbors graph saved.\n')
      } else {
        cat("Warning: RNA_snn graph not found. Skipping saving shared nearest neighbors graph.\n")
      }
      
      
      # Find Clusters
      clustering_resolution = args$clustering_resolution
      clustering_algorithm = args$clustering_algorithm
      
      cat(paste('Starting FindClusters using algorithm:', clustering_algorithm, 'with resolution:', clustering_resolution, '\n'))
      
      algorithm_map = c('louvain' = 1, 'multilevel' = 2, 'slm' = 3, 'leiden' = 4)
      selected_algorithm = algorithm_map[clustering_algorithm]
      
      if (is.null(selected_algorithm) || is.na(selected_algorithm)) {
        cat(paste("Error: Invalid clustering algorithm specified:", clustering_algorithm, ". Please choose from louvain, multilevel, leiden, slm.\n"))
        selected_algorithm = 4 # Default to leiden if invalid
        cat(paste("Defaulting to leiden algorithm (", clustering_algorithm, " -> ", names(algorithm_map)[algorithm_map == selected_algorithm], ").\n"))
      }
      
      seurat_object = FindClusters(object = seurat_object,
                                   resolution = clustering_resolution,
                                   algorithm = selected_algorithm,
                                   verbose = FALSE)
      
      cat("FindClusters complete.\n")
      cat("Current identity class is:", if(length(Idents(seurat_object)) > 0) Idents(seurat_object)[1] else "No cells/identities found", "...\n")
      
      cat('Performing UMAP dimensionality reduction...\n')
      cat("Starting RunUMAP...\n")
      
      # Check if PCA reduction exists and has valid dimensions before running UMAP
      if (!is.null(seurat_object@reductions$pca) && inherits(seurat_object@reductions$pca, "DimReduc") && ncol(seurat_object@reductions$pca@feature.loadings) >= max(valid_pca_dims)) {
        seurat_object = RunUMAP(object = seurat_object,
                                dims = valid_pca_dims, # Use valid_pca_dims
                                verbose = FALSE)
        cat("RunUMAP complete.\n")
        umap_successful = TRUE
      } else {
        cat("Warning: PCA reduction not found or insufficient dimensions for UMAP. Skipping RunUMAP.\n")
        umap_successful = FALSE
      }
      
      
      if (umap_successful) {
        cat('Generating UMAP plot (Static with labels and legend)...\n')
        p9 <- tryCatch({
          DimPlot(object = seurat_object, reduction = 'umap', label = TRUE, label.color = 'black', repel = TRUE, pt.size = 0.7) +
            theme_classic() +
            labs(title = paste0(args$project_name, ' - UMAP Visualization'), subtitle = paste0('Clustering: ', clustering_algorithm, ' (Res: ', clustering_resolution, ') | PCs used: 1-', tail(valid_pca_dims, 1)), x = 'UMAP 1', y = 'UMAP 2') +
            theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 14), plot.subtitle = element_text(hjust = 0.5, size = 10))
        }, error = function(e) {
          cat(paste0("Error generating UMAP plot: ", e$message, "\n"))
          return(NULL)
        })
        
        
        # Save UMAP plot (static) if successfully created
        if (!is.null(p9)) {
          ggsave(filename = paste0(args$project_name, '_UMAP.png'),
                 plot = p9,
                 width = 10,
                 height = 8,
                 units = 'in',
                 dpi = 300)
          cat('UMAP plot (static) saved.\n')
        } else {
          cat('Skipping saving UMAP plot due to previous error.\n')
        }
        
        
      } else {
        cat("Skipping UMAP plot generation due to previous error or warning.\n")
      }
      
    } else {
      cat("Warning: No valid PCA dimensions found. Skipping FindNeighbors, FindClusters, and UMAP.\n")
    }
    
    
    # Find differentially expressed genes based on the argument
    if (args$find_markers) {
      cat('Finding differentially expressed genes...\n')
      cat("Starting FindAllMarkers...\n")
      cat("Using 'DESeq2' for differential gene expression analysis. Note: 'wilcox' or 'negbinom' are often more appropriate for single-cell data.\n")
      
      # Check if clustering identities exist before finding markers
      if (!is.null(Idents(seurat_object)) && length(unique(Idents(seurat_object))) > 1 && ncol(seurat_object) > 0) {
        # Add check for minimum cells per identity for DESeq2 if needed, or switch test.use
        # DESeq2 is sensitive to low cell counts per group. Wilcoxon is more robust.
        # Keeping DESeq2 as requested but adding a note.
        min_cells_per_ident = min(table(Idents(seurat_object)))
        if (min_cells_per_ident < 2) {
          cat(paste0("Warning: Minimum cells per cluster is ", min_cells_per_ident, ". DESeq2 requires at least 2 cells per group. Consider using test.use='wilcox' or test.use='negbinom' instead.\n"))
          # You could add logic here to switch test.use if needed
        }
        
        # Ensure DefaultAssay is RNA for FindAllMarkers
        DefaultAssay(seurat_object) <- "RNA"
        
        markers <- tryCatch({
          FindAllMarkers(object = seurat_object,
                         logfc.threshold = 0.25,
                         min.pct = 0.1,
                         test.use = 'DESeq2',
                         return.thresh = 0.05,
                         verbose = FALSE)
        }, error = function(e) {
          cat(paste0("Error during FindAllMarkers (DESeq2): ", e$message, "\n"))
          cat("Skipping marker finding.\n")
          return(NULL) # Return NULL if FindAllMarkers fails
        })
        
        
        if (!is.null(markers)) {
          write.csv(x = markers,
                    file = paste0(args$project_name, '_markers_res', clustering_resolution, '_alg', clustering_algorithm, '.csv'),
                    row.names = FALSE)
          
          cat('Differential gene expression analysis complete. Markers saved.\n')
        } else {
          cat('Skipping saving markers file due to FindAllMarkers failure.\n')
        }
        
      } else {
        cat("Warning: Clustering identities not found, only one cluster exists, or no cells remain after filtering. Skipping FindAllMarkers.\n")
      }
      
      
    } else {
      cat('Skipping the step of finding differentially expressed genes as per user request.\n')
    }
    
    # Save the Seurat object
    cat('Saving final Seurat object...\n')
    saveRDS(seurat_object, file = paste0(args$project_name, '_processed.rds'))
    cat('Seurat object saving complete.\n')
    
    
    cat('Pipeline finished successfully. Results and plots saved in the working directory.\n')
    
  } else {
    cat("Skipping clustering, UMAP, and marker finding due to previous errors or warnings in scaling or PCA steps.\n")
    cat('Pipeline finished with warnings/errors. Results saved up to that point.\n')
  }
  
  
} else {
  
  cat('Error: Seurat object was not loaded or created properly or has no cells. Pipeline aborted at the start.\n')
  
}
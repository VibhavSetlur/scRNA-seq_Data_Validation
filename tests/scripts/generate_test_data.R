#!/usr/bin/env Rscript

# Test Data Generation Script for snRNA-seq Pipeline
# This script generates synthetic single-cell RNA-seq data for testing

# Set library path to include user-installed packages
.libPaths(c("~/.local/lib/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
  library(dplyr)
  library(stringr)
})

# Function to print colored output
print_status <- function(message, type = "info") {
  colors <- list(
    info = "\033[34m",    # Blue
    success = "\033[32m", # Green
    warning = "\033[33m", # Yellow
    error = "\033[31m",   # Red
    reset = "\033[0m"     # Reset
  )
  cat(paste0(colors[[type]], message, colors$reset, "\n"))
}

# Function to generate synthetic gene names
generate_gene_names <- function(n_genes) {
  # Generate realistic gene names
  gene_prefixes <- c("GENE", "ENSG", "MIR", "LINC", "SNORD", "RP")
  gene_suffixes <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  
  genes <- character(n_genes)
  for (i in 1:n_genes) {
    prefix <- sample(gene_prefixes, 1)
    suffix <- paste0(sample(gene_suffixes, sample(1:3, 1)), collapse = "")
    number <- sample(1000:9999, 1)
    genes[i] <- paste0(prefix, number, suffix)
  }
  
  return(genes)
}

# Function to generate synthetic cell barcodes
generate_cell_barcodes <- function(n_cells) {
  # Generate realistic cell barcodes (similar to 10X Genomics format)
  barcodes <- character(n_cells)
  for (i in 1:n_cells) {
    # Generate 16-character barcode
    barcode <- paste0(sample(c(LETTERS, 0:9), 16, replace = TRUE), collapse = "")
    barcodes[i] <- paste0(barcode, "-1")  # Add -1 suffix like 10X
  }
  
  return(barcodes)
}

# Function to generate synthetic expression matrix with cell types
generate_expression_matrix <- function(n_genes, n_cells, sparsity = 0.9) {
  print_status("Generating synthetic expression matrix with cell types...", "info")
  
  # Set random seed for reproducibility
  set.seed(42)
  
  # Define cell types and their proportions
  cell_types <- c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte")
  cell_proportions <- c(0.4, 0.3, 0.2, 0.1)  # 40% neurons, 30% astrocytes, etc.
  
  # Assign cell types to cells
  cell_type_assignments <- sample(cell_types, n_cells, replace = TRUE, prob = cell_proportions)
  
  # Create marker genes for each cell type
  marker_genes_per_type <- 50  # 50 marker genes per cell type
  total_marker_genes <- length(cell_types) * marker_genes_per_type
  
  # Generate marker gene names
  marker_genes <- character(total_marker_genes)
  for (i in 1:total_marker_genes) {
    cell_type_idx <- ((i - 1) %/% marker_genes_per_type) + 1
    gene_idx <- ((i - 1) %% marker_genes_per_type) + 1
    marker_genes[i] <- paste0(cell_types[cell_type_idx], "_GENE", gene_idx)
  }
  
  # Generate remaining genes (non-marker genes)
  remaining_genes <- n_genes - total_marker_genes
  if (remaining_genes > 0) {
    other_genes <- generate_gene_names(remaining_genes)
    all_gene_names <- c(marker_genes, other_genes)
  } else {
    all_gene_names <- marker_genes
  }
  
  # Create expression matrix with cell type-specific patterns
  matrix_data <- Matrix(0, nrow = n_genes, ncol = n_cells, sparse = TRUE)
  
  # For each cell type, set higher expression for their marker genes
  for (cell_type_idx in 1:length(cell_types)) {
    cell_type <- cell_types[cell_type_idx]
    cells_of_type <- which(cell_type_assignments == cell_type)
    
    if (length(cells_of_type) > 0) {
      # Get marker genes for this cell type
      marker_start <- (cell_type_idx - 1) * marker_genes_per_type + 1
      marker_end <- cell_type_idx * marker_genes_per_type
      cell_type_markers <- marker_start:marker_end
      
      # Set higher expression for marker genes in this cell type
      for (marker_gene_idx in cell_type_markers) {
        if (marker_gene_idx <= n_genes) {
          # Higher expression for marker genes in their cell type
          high_expression <- rnbinom(length(cells_of_type), mu = 20, size = 3)
          matrix_data[marker_gene_idx, cells_of_type] <- high_expression
        }
      }
    }
  }
  
  # Add background expression (sparse)
  total_entries <- n_genes * n_cells
  non_zero_entries <- round(total_entries * (1 - sparsity))
  
  # Generate random positions for additional non-zero entries
  existing_nonzero <- sum(matrix_data > 0)
  additional_entries <- max(0, non_zero_entries - existing_nonzero)
  
  if (additional_entries > 0) {
    # Find positions that are currently zero
    zero_positions <- which(matrix_data == 0)
    if (length(zero_positions) > additional_entries) {
      random_positions <- sample(zero_positions, additional_entries)
      
      # Convert positions to row and column indices
      rows <- ((random_positions - 1) %% n_genes) + 1
      cols <- ((random_positions - 1) %/% n_genes) + 1
      
      # Generate background expression (lower than marker genes)
      background_counts <- rnbinom(additional_entries, mu = 3, size = 2)
      background_counts[background_counts == 0] <- 1
      
      # Add to matrix
      for (i in 1:additional_entries) {
        matrix_data[rows[i], cols[i]] <- background_counts[i]
      }
    }
  }
  
  print_status(paste("Generated matrix with", n_genes, "genes and", n_cells, "cells"), "success")
  print_status(paste("Cell types:", paste(cell_types, collapse = ", ")), "info")
  print_status(paste("Sparsity:", round(1 - (sum(matrix_data > 0) / total_entries), 3)), "info")
  
  # Return matrix and metadata
  return(list(
    matrix = matrix_data,
    gene_names = all_gene_names,
    cell_types = cell_type_assignments,
    marker_genes = marker_genes
  ))
}

# Function to create H5 file
create_h5_file <- function(matrix_data, gene_names, cell_barcodes, output_file) {
  print_status("Creating H5 file...", "info")
  
  # Create H5 file
  h5_file <- H5File$new(output_file, mode = "w")
  
  # Create matrix group
  matrix_group <- h5_file$create_group("matrix")
  
  # Add data
  matrix_group[["data"]] <- matrix_data@x
  matrix_group[["indices"]] <- matrix_data@i
  matrix_group[["indptr"]] <- matrix_data@p
  matrix_group[["shape"]] <- c(nrow(matrix_data), ncol(matrix_data))
  
  # Add feature names (genes)
  features_group <- matrix_group$create_group("features")
  features_group[["name"]] <- gene_names
  
  # Add feature types (all "Gene Expression")
  features_group[["feature_type"]] <- rep("Gene Expression", length(gene_names))
  
  # Add feature genome (all "GRCh38")
  features_group[["genome"]] <- rep("GRCh38", length(gene_names))
  
  # Add barcodes
  matrix_group[["barcodes"]] <- cell_barcodes
  
  # Close file
  h5_file$close()
  
  print_status(paste("H5 file created:", output_file), "success")
}

# Function to create SouporCell test data
create_souporcell_data <- function(cell_barcodes, output_file) {
  print_status("Creating SouporCell test data...", "info")
  
  # Create realistic SouporCell output
  n_cells <- length(cell_barcodes)
  
  # Assign random statuses (mostly singlets, some doublets)
  statuses <- sample(c("singlet", "doublet", "unassigned"), 
                    n_cells, 
                    prob = c(0.85, 0.10, 0.05), 
                    replace = TRUE)
  
  # Assign random clusters
  clusters <- sample(1:5, n_cells, replace = TRUE)
  
  # Create data frame
  souporcell_data <- data.frame(
    barcode = cell_barcodes,
    status = statuses,
    cluster = clusters,
    log_prob_singlet = rnorm(n_cells, mean = -0.5, sd = 1),
    log_prob_doublet = rnorm(n_cells, mean = -2, sd = 1)
  )
  
  # Write to file
  write.table(souporcell_data, 
              file = output_file, 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
  
  print_status(paste("SouporCell data created:", output_file), "success")
}

# Function to create metadata
create_metadata <- function(n_cells, output_file) {
  print_status("Creating metadata...", "info")
  
  # Create realistic metadata
  metadata <- data.frame(
    cell_id = 1:n_cells,
    n_genes = rpois(n_cells, lambda = 1500),
    n_umis = rpois(n_cells, lambda = 5000),
    percent_mt = runif(n_cells, min = 0, max = 15),
    sample_id = rep("test_sample", n_cells),
    batch = rep("batch1", n_cells)
  )
  
  # Write to file
  write.csv(metadata, file = output_file, row.names = FALSE)
  
  print_status(paste("Metadata created:", output_file), "success")
}

# Function to validate generated data
validate_data <- function(matrix_data, gene_names, cell_barcodes) {
  print_status("Validating generated data...", "info")
  
  # Check dimensions
  if (nrow(matrix_data) != length(gene_names)) {
    print_status("Error: Number of genes doesn't match matrix rows", "error")
    return(FALSE)
  }
  
  if (ncol(matrix_data) != length(cell_barcodes)) {
    print_status("Error: Number of cells doesn't match matrix columns", "error")
    return(FALSE)
  }
  
  # Check for valid gene names
  if (any(is.na(gene_names) | gene_names == "")) {
    print_status("Error: Invalid gene names found", "error")
    return(FALSE)
  }
  
  # Check for valid cell barcodes
  if (any(is.na(cell_barcodes) | cell_barcodes == "")) {
    print_status("Error: Invalid cell barcodes found", "error")
    return(FALSE)
  }
  
  # Check matrix properties
  print_status(paste("Matrix dimensions:", nrow(matrix_data), "x", ncol(matrix_data)), "info")
  print_status(paste("Non-zero entries:", length(matrix_data@x)), "info")
  print_status(paste("Sparsity:", round(1 - (length(matrix_data@x) / (nrow(matrix_data) * ncol(matrix_data))), 3)), "info")
  
  print_status("Data validation passed", "success")
  return(TRUE)
}

# Main function
main <- function() {
  print_status("=== Test Data Generation ===", "info")
  
  # Parameters
  n_genes <- 2000
  n_cells <- 1000
  sparsity <- 0.9
  marker_genes_per_type <- 50  # 50 marker genes per cell type
  
  print_status(paste("Generating test data with", n_genes, "genes and", n_cells, "cells"), "info")
  
  # Create output directory
  output_dir <- "tests/data"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate data
  cell_barcodes <- generate_cell_barcodes(n_cells)
  expression_data <- generate_expression_matrix(n_genes, n_cells, sparsity)
  
  # Extract components
  matrix_data <- expression_data$matrix
  gene_names <- expression_data$gene_names
  cell_types <- expression_data$cell_types
  marker_genes <- expression_data$marker_genes
  
  # Validate data
  if (!validate_data(matrix_data, gene_names, cell_barcodes)) {
    print_status("Data validation failed", "error")
    quit(status = 1)
  }
  
  # Create output files
  h5_file <- file.path(output_dir, "test_data.h5")
  souporcell_file <- file.path(output_dir, "test_souporcell.tsv")
  metadata_file <- file.path(output_dir, "test_metadata.csv")
  
  create_h5_file(matrix_data, gene_names, cell_barcodes, h5_file)
  create_souporcell_data(cell_barcodes, souporcell_file)
  create_metadata(n_cells, metadata_file)
  
  # Create summary file
  summary_file <- file.path(output_dir, "test_data_summary.txt")
  cat("Test Data Summary\n", file = summary_file)
  cat("=================\n", file = summary_file, append = TRUE)
  cat("Generated on:", Sys.time(), "\n", file = summary_file, append = TRUE)
  cat("Number of genes:", n_genes, "\n", file = summary_file, append = TRUE)
  cat("Number of cells:", n_cells, "\n", file = summary_file, append = TRUE)
  cat("Matrix sparsity:", round(1 - (sum(matrix_data > 0) / (n_genes * n_cells)), 3), "\n", file = summary_file, append = TRUE)
  cat("Total UMIs:", sum(matrix_data), "\n", file = summary_file, append = TRUE)
  cat("Mean UMIs per cell:", round(mean(colSums(matrix_data)), 1), "\n", file = summary_file, append = TRUE)
  cat("Mean genes per cell:", round(mean(colSums(matrix_data > 0)), 1), "\n", file = summary_file, append = TRUE)
  cat("\nCell Type Distribution:\n", file = summary_file, append = TRUE)
  cell_type_counts <- table(cell_types)
  for (cell_type in names(cell_type_counts)) {
    cat("  ", cell_type, ":", cell_type_counts[cell_type], "cells (", 
        round(cell_type_counts[cell_type]/n_cells*100, 1), "%)\n", 
        file = summary_file, append = TRUE)
  }
  cat("\nMarker Genes:\n", file = summary_file, append = TRUE)
  cat("  Total marker genes:", length(marker_genes), "\n", file = summary_file, append = TRUE)
  cat("  Marker genes per cell type:", marker_genes_per_type, "\n", file = summary_file, append = TRUE)
  
  print_status("=== Test Data Generation Completed ===", "success")
  print_status(paste("Files created in:", output_dir), "info")
  print_status("Files created:", "info")
  print_status(paste("  -", h5_file), "info")
  print_status(paste("  -", souporcell_file), "info")
  print_status(paste("  -", metadata_file), "info")
  print_status(paste("  -", summary_file), "info")
}

# Run main function
main()

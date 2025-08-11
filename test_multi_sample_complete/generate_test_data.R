#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# Function to create mock Seurat object
create_mock_seurat <- function(sample_name, n_cells, n_genes) {
  set.seed(42)
  
  # Create random count matrix
  counts <- matrix(
    rpois(n_cells * n_genes, lambda = 2),
    nrow = n_genes,
    ncol = n_cells
  )
  
  # Create gene and cell names
  gene_names <- paste0("Gene_", 1:n_genes)
  cell_names <- paste0("Cell_", 1:n_cells)
  
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name
  )
  
  # Add metadata
  seurat_obj$sample_id <- sample_name
  seurat_obj$condition <- ifelse(grepl("Control", sample_name), "Control", "Treatment")
  
  return(seurat_obj)
}

# Create test samples
samples <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
n_cells_per_sample <- c(150, 200, 180, 220)
n_genes_per_sample <- c(2500, 3000, 2800, 3200)

for (i in seq_along(samples)) {
  sample_name <- samples[i]
  n_cells <- n_cells_per_sample[i]
  n_genes <- n_genes_per_sample[i]
  
  print(paste("Creating", sample_name, "with", n_cells, "cells and", n_genes, "genes"))
  
  seurat_obj <- create_mock_seurat(sample_name, n_cells, n_genes)
  
  # Save to file
  output_file <- file.path("data", paste0(sample_name, ".rds"))
  saveRDS(seurat_obj, file = output_file)
  
  print(paste("Saved", sample_name, "to", output_file))
}

print("Test data generation completed!")

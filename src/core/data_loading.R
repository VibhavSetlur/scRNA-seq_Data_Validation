# Data Loading Module for snRNA-seq Pipeline
# Handles loading of H5 files and RDS files

# Function to load data from various sources
load_data <- function(args, config) {
  logger_info("Loading data...")
  
  seurat_object <- NULL
  input_type <- NULL
  
  # Load from RDS file
  if (!is.null(args$rds_input) && file.exists(args$rds_input)) {
    logger_info(paste("Loading Seurat object from:", args$rds_input))
    seurat_object <- load_rds_data(args$rds_input)
    input_type <- 'rds'
  }
  
  # Load from H5 file
  else if (!is.null(args$h5_input) && file.exists(args$h5_input)) {
    logger_info(paste("Loading raw data from:", args$h5_input))
    seurat_object <- load_h5_data(args$h5_input, args$project_name)
    input_type <- 'h5'
  }
  
  # Error if no valid input
  else {
    logger_error("No valid input file provided")
    stop("Please provide a valid RDS or H5 input file")
  }
  
  # Validate loaded object
  if (!is.null(seurat_object)) {
    seurat_object <- validate_seurat_object(seurat_object)
    logger_success(paste("Data loaded successfully. Object contains", ncol(seurat_object), "cells and", nrow(seurat_object), "genes"))
  }
  
  return(seurat_object)
}

# Function to load RDS data
load_rds_data <- function(rds_file) {
  tryCatch({
    # Enhanced file validation
    if (!file.exists(rds_file)) {
      stop(paste("RDS file does not exist:", rds_file))
    }
    
    if (file.size(rds_file) == 0) {
      stop(paste("RDS file is empty:", rds_file))
    }
    
    # Check if file is readable
    if (!file.access(rds_file, mode = 4) == 0) {
      stop(paste("RDS file is not readable:", rds_file))
    }
    
    logger_info(paste("Loading RDS file:", rds_file))
    seurat_object <- readRDS(file = rds_file)
    
    # Check if it's a Seurat object
    if (!inherits(seurat_object, "Seurat")) {
      stop("RDS file does not contain a Seurat object")
    }
    
    logger_info("RDS file loaded successfully")
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error loading RDS file:", e$message))
    logger_error(paste("File path:", rds_file))
    logger_error(paste("Error occurred at:", Sys.time()))
    stop(e)
  })
}

# Function to load H5 data
load_h5_data <- function(h5_file, project_name) {
  tryCatch({
    # Enhanced file validation
    if (!file.exists(h5_file)) {
      stop(paste("H5 file does not exist:", h5_file))
    }
    
    if (file.size(h5_file) == 0) {
      stop(paste("H5 file is empty:", h5_file))
    }
    
    # Check if file is readable
    if (!file.access(h5_file, mode = 4) == 0) {
      stop(paste("H5 file is not readable:", h5_file))
    }
    
    # Check if hdf5r package is available
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
      stop("hdf5r package is required to read H5 files. Please install it using: install.packages('hdf5r')")
    }
    
    logger_info(paste("Loading H5 file:", h5_file))
    
    # Read 10X H5 data with error handling
    h5_data <- tryCatch({
      Seurat::Read10X_h5(file = h5_file)
    }, error = function(e) {
      stop(paste("Failed to read H5 file as 10X format:", e$message))
    })
    
    # Validate the loaded data
    if (is.null(h5_data) || length(h5_data) == 0) {
      stop("H5 file contains no valid data")
    }
    
    # Create Seurat object
    seurat_object <- tryCatch({
      Seurat::CreateSeuratObject(
        counts = h5_data,
        project = project_name
      )
    }, error = function(e) {
      stop(paste("Failed to create Seurat object from H5 data:", e$message))
    })
    
    logger_info("H5 file loaded and Seurat object created successfully")
    return(seurat_object)
    
  }, error = function(e) {
    logger_error(paste("Error loading H5 file:", e$message))
    logger_error(paste("File path:", h5_file))
    logger_error(paste("Error occurred at:", Sys.time()))
    stop(e)
  })
}

# Function to validate Seurat object
validate_seurat_object <- function(seurat_object) {
  logger_info("Validating Seurat object...")
  
  # Check basic structure
  if (!inherits(seurat_object, "Seurat")) {
    stop("Object is not a Seurat object")
  }
  
  # Check for required slots
  required_slots <- c("assays", "meta.data", "active.assay")
  missing_slots <- setdiff(required_slots, names(seurat_object@.Data))
  
  if (length(missing_slots) > 0) {
    logger_warning(paste("Missing slots in Seurat object:", paste(missing_slots, collapse = ", ")))
  }
  
  # Check for RNA assay
  if (!"RNA" %in% names(seurat_object@assays)) {
    logger_warning("RNA assay not found, creating default assay")
    if (length(names(seurat_object@assays)) > 0) {
      Seurat::DefaultAssay(seurat_object) <- names(seurat_object@assays)[1]
    }
  } else {
    Seurat::DefaultAssay(seurat_object) <- "RNA"
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
  
  logger_info(paste("Validation passed. Object dimensions:", n_genes, "genes x", n_cells, "cells"))
  
  return(seurat_object)
}

# Function to load SouporCell data
load_souporcell_data <- function(souporcell_file) {
  if (is.null(souporcell_file) || !file.exists(souporcell_file)) {
    logger_warning("SouporCell file not provided or not found")
    return(NULL)
  }
  
  tryCatch({
    logger_info(paste("Loading SouporCell data from:", souporcell_file))
    
    # Read SouporCell output
    souporcell_data <- read.csv(
      souporcell_file,
      sep = '\t',
      stringsAsFactors = FALSE
    )
    
    # Validate required columns
    required_columns <- c("barcode", "status")
    missing_columns <- setdiff(required_columns, colnames(souporcell_data))
    
    if (length(missing_columns) > 0) {
      stop(paste("SouporCell file missing required columns:", paste(missing_columns, collapse = ", ")))
    }
    
    logger_success(paste("SouporCell data loaded successfully. Found", nrow(souporcell_data), "cells"))
    
    return(souporcell_data)
    
  }, error = function(e) {
    logger_error(paste("Error loading SouporCell data:", e$message))
    stop(e)
  })
}

# Function to get data summary
get_data_summary <- function(seurat_object) {
  summary <- list(
    n_cells = ncol(seurat_object),
    n_genes = nrow(seurat_object),
    n_features = sum(seurat_object$nFeature_RNA),
    n_counts = sum(seurat_object$nCount_RNA),
    mean_features_per_cell = mean(seurat_object$nFeature_RNA),
    mean_counts_per_cell = mean(seurat_object$nCount_RNA),
    median_features_per_cell = median(seurat_object$nFeature_RNA),
    median_counts_per_cell = median(seurat_object$nCount_RNA)
  )
  
  return(summary)
}

# Function to print data summary
print_data_summary <- function(seurat_object) {
  summary <- get_data_summary(seurat_object)
  
  logger_info("Data Summary:")
  logger_info(paste("  Number of cells:", summary$n_cells))
  logger_info(paste("  Number of genes:", summary$n_genes))
  logger_info(paste("  Total features:", summary$n_features))
  logger_info(paste("  Total counts:", summary$n_counts))
  logger_info(paste("  Mean features per cell:", round(summary$mean_features_per_cell, 1)))
  logger_info(paste("  Mean counts per cell:", round(summary$mean_counts_per_cell, 1)))
  logger_info(paste("  Median features per cell:", round(summary$median_features_per_cell, 1)))
  logger_info(paste("  Median counts per cell:", round(summary$median_counts_per_cell, 1)))
}

# Function to save data summary
save_data_summary <- function(seurat_object, output_dir, project_name) {
  summary <- get_data_summary(seurat_object)
  
  # Create summary data frame
  summary_df <- data.frame(
    Metric = names(summary),
    Value = unlist(summary)
  )
  
  # Save to file
  summary_file <- file.path(output_dir, paste0(project_name, '_data_summary.csv'))
  write.csv(summary_df, file = summary_file, row.names = FALSE)
  
  logger_info(paste("Data summary saved to:", summary_file))
  
  return(summary_file)
}

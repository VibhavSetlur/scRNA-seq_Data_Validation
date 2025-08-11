#!/usr/bin/env Rscript

# snRNA-seq Pipeline Setup Script
# This script installs all required dependencies and checks system requirements

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    install.packages("yaml")
  }
  library(yaml)
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

# Function to check R version
check_r_version <- function() {
  print_status("Checking R version...", "info")
  r_version <- R.version.string
  print_status(paste("R version:", r_version), "info")
  
  # Check if R version is >= 4.3.0
  major_version <- as.numeric(R.version$major)
  minor_version <- as.numeric(R.version$minor)
  
  if (major_version < 4 || (major_version == 4 && minor_version < 3)) {
    print_status("Warning: R version 4.3.0 or higher is recommended", "warning")
  } else {
    print_status("✓ R version is compatible", "success")
  }
}

# Function to check system resources
check_system_resources <- function() {
  print_status("Checking system resources...", "info")
  
  # Check available memory
  if (Sys.info()["sysname"] == "Linux") {
    mem_info <- system("free -h", intern = TRUE)
    print_status("Memory information:", "info")
    cat(paste(mem_info, collapse = "\n"), "\n")
  } else if (Sys.info()["sysname"] == "Darwin") {
    mem_info <- system("vm_stat", intern = TRUE)
    print_status("Memory information (macOS):", "info")
    cat(paste(mem_info[1:10], collapse = "\n"), "\n")
  }
  
  # Check available disk space
  disk_space <- file.size(".") / (1024^3)  # GB
  print_status(paste("Available disk space in current directory:", round(disk_space, 2), "GB"), "info")
  
  # Check number of CPU cores
  n_cores <- parallel::detectCores()
  print_status(paste("Number of CPU cores:", n_cores), "info")
}

# Function to install packages from requirements file
install_packages_from_requirements <- function() {
  print_status("Installing packages from requirements.txt...", "info")
  
  # Read requirements file
  if (!file.exists("requirements.txt")) {
    print_status("Error: requirements.txt not found", "error")
    return(FALSE)
  }
  
  requirements <- readLines("requirements.txt")
  requirements <- requirements[!grepl("^#", requirements) & requirements != ""]
  
  # Extract package names and versions
  packages <- character()
  for (req in requirements) {
    if (grepl(">=", req)) {
      pkg_name <- trimws(strsplit(req, ">=")[[1]][1])
      packages <- c(packages, pkg_name)
    } else {
      pkg_name <- trimws(req)
      packages <- c(packages, pkg_name)
    }
  }
  
  # Install packages only if not already installed
  for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      print_status(paste("✓", pkg, "is already installed"), "success")
    } else {
      print_status(paste("Installing", pkg, "..."), "info")
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        if (requireNamespace(pkg, quietly = TRUE)) {
          print_status(paste("✓ Successfully installed", pkg), "success")
        } else {
          print_status(paste("Failed to install", pkg), "error")
        }
      }, error = function(e) {
        print_status(paste("Error installing", pkg, ":", e$message), "error")
      })
    }
  }
  
  return(TRUE)
}

# Function to install Bioconductor packages
install_bioconductor_packages <- function() {
  print_status("Installing Bioconductor packages...", "info")
  
  # Check if BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    print_status("Installing BiocManager...", "info")
    install.packages("BiocManager")
  } else {
    print_status("✓ BiocManager is already installed", "success")
  }
  
  bioc_packages <- c("Seurat", "DESeq2", "edgeR", "sparseMatrixStats")
  
  for (pkg in bioc_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      print_status(paste("✓", pkg, "is already installed"), "success")
    } else {
      print_status(paste("Installing Bioconductor package:", pkg), "info")
      tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        if (requireNamespace(pkg, quietly = TRUE)) {
          print_status(paste("✓ Successfully installed", pkg), "success")
        } else {
          print_status(paste("Failed to install", pkg), "error")
        }
      }, error = function(e) {
        print_status(paste("Error installing", pkg, ":", e$message), "error")
      })
    }
  }
}

# Function to install additional packages for multi-sample processing and Shiny
install_additional_packages <- function() {
  print_status("Installing additional packages for multi-sample processing and Shiny...", "info")
  
  additional_packages <- c(
    "future", "future.apply", "parallel", "purrr", "shinyFiles",
    "shinyWidgets", "shinyjs", "DT", "fs", "promises", "Matrix", 
    "hdf5r", "stringr", "markdown"
  )
  
  for (pkg in additional_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      print_status(paste("✓", pkg, "is already installed"), "success")
    } else {
      print_status(paste("Installing", pkg, "..."), "info")
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        if (requireNamespace(pkg, quietly = TRUE)) {
          print_status(paste("✓ Successfully installed", pkg), "success")
        } else {
          print_status(paste("Failed to install", pkg), "error")
        }
      }, error = function(e) {
        print_status(paste("Error installing", pkg, ":", e$message), "error")
      })
    }
  }
}

# Function to create necessary directories
create_directories <- function() {
  print_status("Creating necessary directories...", "info")
  
  dirs <- c(
    "outputs",
    "logs", 
    "temp",
    "data",
    "results"
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      print_status(paste("✓ Created directory:", dir), "success")
    } else {
      print_status(paste("✓ Directory already exists:", dir), "success")
    }
  }
}

# Function to validate installation
validate_installation <- function() {
  print_status("Validating installation...", "info")
  
  required_packages <- c(
    "Seurat", "tidyverse", "patchwork", "argparse", "scales", "ggrepel",
    "shiny", "shinydashboard", "plotly", "yaml", "testthat", "shinyFiles",
    "shinyWidgets", "shinyjs", "DT", "fs", "future", "promises", "dplyr",
    "Matrix", "hdf5r", "stringr"
  )
  
  missing_packages <- character()
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) == 0) {
    print_status("✓ All required packages are installed", "success")
    return(TRUE)
  } else {
    print_status("Missing packages:", "error")
    for (pkg in missing_packages) {
      print_status(paste("  -", pkg), "error")
    }
    return(FALSE)
  }
}

# Function to show more information
show_more_info <- function() {
  print_status("\n=== SHOW MORE INFORMATION ===", "info")
  
  print_status("System Information:", "info")
  cat("Operating System:", Sys.info()["sysname"], Sys.info()["release"], "\n")
  cat("Machine:", Sys.info()["machine"], "\n")
  cat("User:", Sys.info()["user"], "\n")
  cat("Working Directory:", getwd(), "\n")
  
  print_status("\nR Session Information:", "info")
  cat("R Version:", R.version.string, "\n")
  cat("Platform:", R.version$platform, "\n")
  cat("Architecture:", R.version$arch, "\n")
  
  print_status("\nInstalled Packages:", "info")
  installed_packages <- installed.packages()
  cat("Total packages:", nrow(installed_packages), "\n")
  
  print_status("\nBioconductor Version:", "info")
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    cat("BiocManager version:", as.character(packageVersion("BiocManager")), "\n")
  } else {
    cat("BiocManager not installed\n")
  }
  
  print_status("\nSeurat Version:", "info")
  if (requireNamespace("Seurat", quietly = TRUE)) {
    cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
  } else {
    cat("Seurat not installed\n")
  }
  
  print_status("\nMemory Usage:", "info")
  mem_usage <- gc()
  cat("Memory used:", round(sum(mem_usage[,2]) / 1024^2, 2), "MB\n")
  
  print_status("\nAvailable CRAN Mirrors:", "info")
  mirrors <- getCRANmirrors()
  cat("Number of mirrors:", nrow(mirrors), "\n")
  
  print_status("\nR Library Paths:", "info")
  cat("Library paths:\n")
  for (path in .libPaths()) {
    cat("  -", path, "\n")
  }
}

# Main setup function
main <- function() {
  print_status("=== snRNA-seq Pipeline Setup ===", "info")
  print_status("Starting installation and system check...", "info")
  
  # Check R version
  check_r_version()
  
  # Check system resources
  check_system_resources()
  
  # Install CRAN packages
  install_packages_from_requirements()
  
  # Install Bioconductor packages
  install_bioconductor_packages()
  
  # Install additional packages for multi-sample processing and Shiny
  install_additional_packages()
  
  # Create directories
  create_directories()
  
  # Validate installation
  success <- validate_installation()
  
  if (success) {
    print_status("\n=== SETUP COMPLETED SUCCESSFULLY ===", "success")
    print_status("You can now run the snRNA-seq pipeline!", "success")
  } else {
    print_status("\n=== SETUP COMPLETED WITH ERRORS ===", "error")
    print_status("Please check the error messages above and try again.", "error")
  }
  
  # Show more information
  show_more_info()
  
  print_status("\nSetup script completed!", "info")
}

# Run main function
main()

# Logger utility for snRNA-seq Pipeline
# Provides consistent logging functionality across the pipeline

# Global logger state
logger_state <- list(
  verbose = TRUE,
  log_file = NULL,
  start_time = NULL
)

# Function to initialize logger
init_logger <- function(verbose = TRUE, log_file = NULL) {
  logger_state$verbose <<- verbose
  logger_state$start_time <<- Sys.time()
  
  if (!is.null(log_file)) {
    logger_state$log_file <<- log_file
    # Create log directory if it doesn't exist
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE)
    }
  }
}

# Function to get timestamp
get_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

# Function to get elapsed time
get_elapsed_time <- function() {
  if (is.null(logger_state$start_time)) {
    return("0s")
  }
  
  elapsed <- difftime(Sys.time(), logger_state$start_time, units = "secs")
  if (elapsed < 60) {
    return(paste0(round(elapsed, 1), "s"))
  } else if (elapsed < 3600) {
    return(paste0(round(elapsed / 60, 1), "m"))
  } else {
    return(paste0(round(elapsed / 3600, 1), "h"))
  }
}

# Function to write log message
write_log <- function(level, message, color = NULL) {
  timestamp <- get_timestamp()
  elapsed <- get_elapsed_time()
  log_entry <- paste0("[", timestamp, "] [", elapsed, "] [", level, "] ", message)
  
  # Write to console
  if (logger_state$verbose) {
    cat(log_entry, "\n")
  }
  
  # Write to log file if specified
  if (!is.null(logger_state$log_file)) {
    write(log_entry, logger_state$log_file, append = TRUE)
  }
}

# Logging functions
logger_info <- function(message) {
  write_log("INFO", message, "blue")
}

logger_success <- function(message) {
  write_log("SUCCESS", message, "green")
}

logger_warning <- function(message) {
  write_log("WARNING", message, "yellow")
}

logger_error <- function(message) {
  write_log("ERROR", message, "red")
}

logger_debug <- function(message) {
  if (logger_state$verbose) {
    write_log("DEBUG", message, "magenta")
  }
}

# Function to log memory usage
logger_memory <- function() {
  if (requireNamespace("pryr", quietly = TRUE)) {
    mem_usage <- pryr::mem_used()
    logger_debug(paste("Memory usage:", format(mem_usage, units = "MB")))
  }
}

# Function to log system information
logger_system_info <- function() {
  logger_info("System Information:")
  logger_info(paste("  OS:", Sys.info()["sysname"], Sys.info()["release"]))
  logger_info(paste("  Machine:", Sys.info()["machine"]))
  logger_info(paste("  User:", Sys.info()["user"]))
  logger_info(paste("  R Version:", R.version.string))
  logger_info(paste("  Working Directory:", getwd()))
  
  # Memory information
  if (Sys.info()["sysname"] == "Linux") {
    mem_info <- system("free -h | grep Mem", intern = TRUE)
    if (length(mem_info) > 0) {
      logger_info(paste("  Memory:", mem_info[1]))
    }
  }
  
  # CPU cores
  n_cores <- parallel::detectCores()
  logger_info(paste("  CPU Cores:", n_cores))
}

# Function to create progress bar
create_progress_bar <- function(total, title = "Progress") {
  if (logger_state$verbose) {
    logger_info(paste("Starting:", title))
  }
  
  list(
    total = total,
    current = 0,
    title = title,
    update = function(increment = 1) {
      .self <- parent.env(environment())
      .self$current <- .self$current + increment
      if (logger_state$verbose) {
        percent <- round((.self$current / .self$total) * 100, 1)
        cat(sprintf("\r%s: %d/%d (%s%%)", 
                   .self$title, .self$current, .self$total, percent))
        if (.self$current >= .self$total) {
          cat("\n")
          logger_success(paste("Completed:", .self$title))
        }
      }
    }
  )
}

# Function to log function entry/exit
log_function <- function(func_name) {
  logger_debug(paste("Entering function:", func_name))
  
  function(...) {
    args <- list(...)
    logger_debug(paste("Function", func_name, "called with", length(args), "arguments"))
    
    result <- tryCatch({
      do.call(func_name, args)
    }, error = function(e) {
      logger_error(paste("Error in", func_name, ":", e$message))
      stop(e)
    })
    
    logger_debug(paste("Function", func_name, "completed successfully"))
    return(result)
  }
}

# Function to create summary log
create_summary_log <- function(output_dir, project_name) {
  summary_file <- file.path(output_dir, paste0(project_name, '_pipeline_summary.log'))
  
  sink(summary_file)
  cat("snRNA-seq Pipeline Summary\n")
  cat("==========================\n")
  cat("Project:", project_name, "\n")
  cat("Start Time:", format(logger_state$start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Duration:", get_elapsed_time(), "\n")
  cat("Working Directory:", getwd(), "\n")
  cat("Output Directory:", output_dir, "\n")
  cat("\n")
  
  # System information
  cat("System Information:\n")
  cat("  OS:", Sys.info()["sysname"], Sys.info()["release"], "\n")
  cat("  R Version:", R.version.string, "\n")
  cat("  Platform:", R.version$platform, "\n")
  
  sink()
  
  logger_info(paste("Pipeline summary saved to:", summary_file))
}

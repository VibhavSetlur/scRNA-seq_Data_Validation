#!/usr/bin/env Rscript

# Run Shiny App for snRNA-seq Pipeline
cat("=== Starting snRNA-seq Pipeline Shiny App ===\n")

# Set library path
.libPaths(c("~/.local/lib/R/library", .libPaths()))

# Load required packages
cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(shinyjs)
  library(shinyFiles)
  library(DT)
  library(plotly)
  library(yaml)
  library(fs)
  library(future)
  library(promises)
})

cat("✓ All packages loaded successfully!\n")

# Source the app.R file to get ui and server functions
cat("Loading Shiny application...\n")
source("shiny_app/app.R")

# Get host and port from environment variables or use defaults
host <- Sys.getenv("SHINY_HOST", "0.0.0.0")  # Use 0.0.0.0 for server deployment
port <- as.numeric(Sys.getenv("SHINY_PORT", "3838"))

# Check if port is available, if not try alternative ports
check_port <- function(port) {
  tryCatch({
    con <- socketConnection(host = "localhost", port = port, server = TRUE)
    close(con)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

if (!check_port(port)) {
  cat("Warning: Port", port, "is already in use.\n")
  # Try alternative ports
  alternative_ports <- c(3839, 3840, 3841, 8080, 8081)
  for (alt_port in alternative_ports) {
    if (check_port(alt_port)) {
      port <- alt_port
      cat("Using alternative port:", port, "\n")
      break
    }
  }
}

cat("Starting Shiny app...\n")
cat("Host:", host, "\n")
cat("Port:", port, "\n")
cat("Access the app at: http://", host, ":", port, "\n", sep="")
cat("Press Ctrl+C to stop the app.\n\n")

# Run the Shiny app with server configuration
shinyApp(
  ui = ui, 
  server = server,
  options = list(
    host = host,
    port = port,
    launch.browser = FALSE  # Don't launch browser automatically on server
  )
)

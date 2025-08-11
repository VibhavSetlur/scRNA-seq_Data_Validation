#!/usr/bin/env Rscript

# Test script for Shiny app
cat("=== Testing Shiny App Dependencies ===\n")

# Set library path
.libPaths(c("~/.local/lib/R/library", .libPaths()))

# Test loading packages
cat("Loading packages...\n")
tryCatch({
  library(shiny)
  cat("✓ Shiny loaded\n")
}, error = function(e) {
  cat("✗ Error loading shiny:", e$message, "\n")
})

tryCatch({
  library(shinydashboard)
  cat("✓ Shinydashboard loaded\n")
}, error = function(e) {
  cat("✗ Error loading shinydashboard:", e$message, "\n")
})

tryCatch({
  library(shinyWidgets)
  cat("✓ ShinyWidgets loaded\n")
}, error = function(e) {
  cat("✗ Error loading shinyWidgets:", e$message, "\n")
})

tryCatch({
  library(DT)
  cat("✓ DT loaded\n")
}, error = function(e) {
  cat("✗ Error loading DT:", e$message, "\n")
})

tryCatch({
  library(plotly)
  cat("✓ Plotly loaded\n")
}, error = function(e) {
  cat("✗ Error loading plotly:", e$message, "\n")
})

cat("\n=== Testing Shiny App Script ===\n")
cat("Checking if app.R can be sourced...\n")

tryCatch({
  # Change to shiny_app directory
  setwd("shiny_app")
  source("app.R")
  cat("✓ Shiny app script loaded successfully!\n")
}, error = function(e) {
  cat("✗ Error loading app.R:", e$message, "\n")
  cat("Error details:", e$call, "\n")
})

cat("\n=== Test Complete ===\n")

#!/usr/bin/env Rscript

# snRNA-seq Pipeline Shiny Application
# Provides a user-friendly web interface for the pipeline

# Enable comprehensive debugging and error tracking
options(
  shiny.error = browser,
  shiny.reactlog = TRUE,
  shiny.trace = TRUE,
  warn = 1
)

# Add terminal logging function
log_to_terminal <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, message))
}

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
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# Source utility functions
source("src/utils/logger.R")
source("src/utils/config.R")
source("src/core/pipeline.R")

# UI Definition
ui <- dashboardPage(
  dashboardHeader(
    title = "snRNA-seq Pipeline",
    titleWidth = 300
  ),
  
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Data Input", tabName = "data_input", icon = icon("upload")),
      menuItem("Quality Control", tabName = "qc", icon = icon("check-circle")),
      menuItem("Processing", tabName = "processing", icon = icon("cogs")),
      menuItem("Clustering", tabName = "clustering", icon = icon("sitemap")),
      menuItem("Visualization", tabName = "visualization", icon = icon("chart-bar")),
      menuItem("Run Pipeline", tabName = "run_pipeline", icon = icon("play")),
      menuItem("Results", tabName = "results", icon = icon("folder-open")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "css/style.css"),
      tags$script(src = "js/app.js")
    ),
    
    tabItems(
      # Dashboard Tab
      tabItem(
        tabName = "dashboard",
        fluidRow(
          column(12,
            box(
              title = "Welcome to snRNA-seq Pipeline",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              p("This application provides a user-friendly interface for single-nucleus RNA sequencing analysis using Seurat."),
              p("Follow the tabs on the left to configure and run your analysis.")
            )
          )
        ),
        fluidRow(
          column(6,
            box(
              title = "System Status",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              verbatimTextOutput("system_status")
            )
          ),
          column(6,
            box(
              title = "Recent Projects",
              width = 12,
              status = "success",
              solidHeader = TRUE,
              tableOutput("recent_projects")
            )
          )
        )
      ),
      
      # Data Input Tab
      tabItem(
        tabName = "data_input",
        fluidRow(
          column(12,
            box(
              title = "Multi-Sample Data Input Configuration",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              
              # Input type selection
              radioButtons("input_type", "Input Type:",
                          choices = c("H5 Files (10X Genomics)" = "h5",
                                    "RDS Files (Seurat Objects)" = "rds"),
                          selected = "h5"),
              
              # Multiple file input
              conditionalPanel(
                condition = "input.input_type == 'h5'",
                fileInput("h5_files", "Select H5 Files (Multiple):",
                         accept = c(".h5", ".hdf5"),
                         multiple = TRUE)
              ),
              
              conditionalPanel(
                condition = "input.input_type == 'rds'",
                fileInput("rds_files", "Select RDS Files (Multiple):",
                         accept = ".rds",
                         multiple = TRUE)
              ),
              
              # Project settings
              textInput("project_name", "Project Name:", value = "MultiSampleProject"),
              textInput("working_dir", "Working Directory:", value = "."),
              
              # Parallel processing settings
              checkboxInput("use_parallel", "Use Parallel Processing", value = TRUE),
              conditionalPanel(
                condition = "input.use_parallel == true",
                numericInput("n_cores", "Number of Cores:",
                            value = 2, min = 1, max = 4)
              ),
              
              # SouporCell input
              checkboxInput("use_souporcell", "Use SouporCell for doublet detection"),
              conditionalPanel(
                condition = "input.use_souporcell == true",
                fileInput("souporcell_files", "SouporCell Output Files (Multiple):",
                         accept = c(".tsv", ".txt"),
                         multiple = TRUE)
              ),
              
              # Sample configuration
              checkboxInput("use_sample_configs", "Use Sample-Specific Configurations"),
              conditionalPanel(
                condition = "input.use_sample_configs == true",
                fileInput("sample_configs_file", "Sample Configuration File:",
                         accept = c(".yaml", ".yml"))
              )
            )
          )
        ),
        
        # Sample preview section
        fluidRow(
          column(12,
            box(
              title = "Sample Preview",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              
              verbatimTextOutput("sample_preview"),
              
              actionButton("refresh_samples", "Refresh Sample List", 
                          class = "btn-info")
            )
          )
        )
      ),
      
      # Quality Control Tab
      tabItem(
        tabName = "qc",
        fluidRow(
          column(12,
            box(
              title = "Quality Control Parameters",
              width = 12,
              status = "warning",
              solidHeader = TRUE,
              
              fluidRow(
                column(6,
                  numericInput("min_features", "Minimum Features per Cell:",
                              value = 200, min = 0, max = 10000),
                  numericInput("min_counts", "Minimum Counts per Cell:",
                              value = 1000, min = 0, max = 100000),
                  numericInput("max_features", "Maximum Features per Cell:",
                              value = 6000, min = 0, max = 50000)
                ),
                column(6,
                  numericInput("max_counts", "Maximum Counts per Cell:",
                              value = 25000, min = 0, max = 1000000),
                  numericInput("max_mt_percent", "Maximum Mitochondrial %:",
                              value = 20, min = 0, max = 100),
                  numericInput("min_cells", "Minimum Cells per Gene:",
                              value = 3, min = 0, max = 100)
                )
              )
            )
          )
        )
      ),
      
      # Processing Tab
      tabItem(
        tabName = "processing",
        fluidRow(
          column(12,
            box(
              title = "Data Processing Parameters",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              
              fluidRow(
                column(6,
                  selectInput("normalization_method", "Normalization Method:",
                             choices = c("CLR" = "CLR",
                                       "LogNormalize" = "LogNormalize",
                                       "RC" = "RC"),
                             selected = "CLR"),
                  numericInput("n_variable_features", "Number of Variable Features:",
                              value = 2000, min = 100, max = 10000),
                  selectInput("scaling_method", "Scaling Method:",
                             choices = c("negbinom" = "negbinom",
                                       "linear" = "linear"),
                             selected = "negbinom")
                ),
                column(6,
                  numericInput("pca_dimensions", "PCA Dimensions:",
                              value = 15, min = 1, max = 100),
                  numericInput("scale_factor", "Scale Factor:",
                              value = 10000, min = 1000, max = 100000),
                  checkboxInput("save_intermediate", "Save Intermediate Results",
                               value = TRUE)
                )
              )
            )
          )
        )
      ),
      
      # Clustering Tab
      tabItem(
        tabName = "clustering",
        fluidRow(
          column(12,
            box(
              title = "Clustering Parameters",
              width = 12,
              status = "success",
              solidHeader = TRUE,
              
              fluidRow(
                column(6,
                  selectInput("clustering_algorithm", "Clustering Algorithm:",
                             choices = c("leiden" = "leiden",
                                       "louvain" = "louvain",
                                       "multilevel" = "multilevel",
                                       "slm" = "slm"),
                             selected = "leiden"),
                  numericInput("clustering_resolution", "Clustering Resolution:",
                              value = 0.5, min = 0.1, max = 2.0, step = 0.1),
                  numericInput("min_cluster_size", "Minimum Cluster Size:",
                              value = 10, min = 1, max = 1000)
                ),
                column(6,
                  numericInput("umap_n_neighbors", "UMAP n_neighbors:",
                              value = 30, min = 5, max = 100),
                  numericInput("umap_min_dist", "UMAP min_dist:",
                              value = 0.3, min = 0.01, max = 1.0, step = 0.01),
                  numericInput("random_seed", "Random Seed:",
                              value = 42, min = 1, max = 10000)
                )
              )
            )
          )
        )
      ),
      
      # Visualization Tab
      tabItem(
        tabName = "visualization",
        fluidRow(
          column(12,
            box(
              title = "Visualization Settings",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              
              fluidRow(
                column(6,
                  selectInput("plot_theme", "Plot Theme:",
                             choices = c("classic" = "classic",
                                       "minimal" = "minimal",
                                       "bw" = "bw"),
                             selected = "classic"),
                  selectInput("color_palette", "Color Palette:",
                             choices = c("viridis" = "viridis",
                                       "rainbow" = "rainbow",
                                       "brewer" = "brewer"),
                             selected = "viridis"),
                  numericInput("point_size", "Point Size:",
                              value = 0.7, min = 0.1, max = 5.0, step = 0.1)
                ),
                column(6,
                  numericInput("label_size", "Label Size:",
                              value = 3, min = 1, max = 10),
                  numericInput("title_size", "Title Size:",
                              value = 14, min = 8, max = 24),
                  selectInput("output_format", "Output Format:",
                             choices = c("png" = "png",
                                       "pdf" = "pdf",
                                       "svg" = "svg"),
                             selected = "png")
                )
              )
            )
          )
        )
      ),
      
      # Run Pipeline Tab
      tabItem(
        tabName = "run_pipeline",
        fluidRow(
          column(12,
            box(
              title = "Pipeline Execution",
              width = 12,
              status = "danger",
              solidHeader = TRUE,
              
              # Analysis options
              checkboxInput("find_markers", "Find Cluster Markers", value = TRUE),
              checkboxInput("verbose", "Verbose Output", value = TRUE),
              
              # Run button
              actionButton("run_pipeline", "Run Pipeline", 
                          class = "btn-primary btn-lg",
                          style = "width: 200px; height: 50px; font-size: 18px;"),
              
              # Status indicator
              conditionalPanel(
                condition = "input.run_pipeline > 0",
                div(
                  style = "margin: 10px 0;",
                  uiOutput("pipeline_status_indicator")
                )
              ),
              
              # Progress
              conditionalPanel(
                condition = "input.run_pipeline > 0",
                progressBar(id = "pipeline_progress", value = 0, display_pct = TRUE),
                verbatimTextOutput("pipeline_log")
              )
            )
          )
        )
      ),
      
      # Results Tab
      tabItem(
        tabName = "results",
        fluidRow(
          column(12,
            box(
              title = "Multi-Sample Analysis Results",
              width = 12,
              status = "success",
              solidHeader = TRUE,
              
              tabsetPanel(
                tabPanel("Project Summary", 
                  verbatimTextOutput("project_summary")
                ),
                tabPanel("Sample Selection",
                  fluidRow(
                    column(6,
                      selectInput("selected_sample", "Select Sample:",
                                 choices = NULL)
                    ),
                    column(6,
                      actionButton("load_sample_results", "Load Sample Results",
                                  class = "btn-primary")
                    )
                  ),
                  verbatimTextOutput("sample_summary")
                ),
                tabPanel("Sample Comparisons",
                  uiOutput("comparison_plots")
                ),
                tabPanel("Individual Sample Results",
                  tabsetPanel(
                    tabPanel("Plots",
                      uiOutput("sample_plot_outputs")
                    ),
                    tabPanel("Tables",
                      DT::dataTableOutput("sample_results_table")
                    ),
                    tabPanel("Files",
                      uiOutput("sample_file_outputs")
                    )
                  )
                ),
                tabPanel("Project Files",
                  uiOutput("project_file_outputs")
                )
              )
            )
          )
        )
      ),
      
      # Help Tab
      tabItem(
        tabName = "help",
        fluidRow(
          column(12,
            box(
              title = "Help and Documentation",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              
              tabsetPanel(
                tabPanel("Getting Started",
                  includeMarkdown("shiny_app/docs/getting_started.md")
                ),
                tabPanel("Parameters",
                  includeMarkdown("shiny_app/docs/parameters.md")
                ),
                tabPanel("Troubleshooting",
                  includeMarkdown("shiny_app/docs/troubleshooting.md")
                ),
                tabPanel("FAQ",
                  includeMarkdown("shiny_app/docs/faq.md")
                )
              )
            )
          )
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    pipeline_running = FALSE,
    results = NULL,
    log_messages = character(),
    progress = 0,
    status = "idle",  # idle, running, completed, error
    current_sample_seurat = NULL,
    current_sample_name = NULL
  )
  
  # System status
  output$system_status <- renderPrint({
    cat("R Version:", R.version.string, "\n")
    cat("Platform:", R.version$platform, "\n")
    cat("Working Directory:", getwd(), "\n")
    cat("Available Memory:", format(memory.size(), units = "MB"), "\n")
    cat("CPU Cores:", parallel::detectCores(), "\n")
  })
  
  # Recent projects - Remove dummy data and show actual project information
  output$recent_projects <- renderTable({
    # Log to terminal for debugging
    log_to_terminal("=== Loading Recent Projects ===", "INFO")
    
    # Check for actual project outputs in multiple directories
    project_dirs <- character()
    
    # Check results directory
    results_dir <- "results"
    log_to_terminal(paste("Checking results directory:", results_dir), "INFO")
    if (dir.exists(results_dir)) {
      log_to_terminal("Results directory exists", "INFO")
      subdirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
      subdirs <- subdirs[subdirs != ""]
      log_to_terminal(paste("Found subdirectories in results:", paste(subdirs, collapse = ", ")), "INFO")
    }
    
    # Check for TestProject_outputs (actual pipeline run)
    test_project_dir <- "TestProject_outputs"
    log_to_terminal(paste("Checking test project directory:", test_project_dir), "INFO")
    if (dir.exists(test_project_dir)) {
      log_to_terminal("TestProject_outputs directory exists", "INFO")
      project_dirs <- c(project_dirs, test_project_dir)
    }
    
    # Check for any other *_outputs directories (pipeline outputs)
    all_dirs <- list.dirs(".", full.names = FALSE, recursive = FALSE)
    output_dirs <- all_dirs[grepl("_outputs$", all_dirs)]
    log_to_terminal(paste("Found output directories:", paste(output_dirs, collapse = ", ")), "INFO")
    project_dirs <- c(project_dirs, output_dirs)
    
    # Remove duplicates
    project_dirs <- unique(project_dirs)
    
    log_to_terminal(paste("Total project directories found:", length(project_dirs)), "INFO")
    
    if (length(project_dirs) > 0) {
      log_to_terminal("Processing project directories:", "INFO")
      for (project_dir in project_dirs) {
        log_to_terminal(paste("  -", project_dir), "INFO")
      }
      
      # Get project information
      project_info <- lapply(project_dirs, function(project_dir) {
        # Extract project name from directory name
        project_name <- gsub("_outputs$", "", project_dir)
        
        # Check for session info file to determine if project is completed
        session_file <- file.path(project_dir, paste0(project_name, "_session_info.txt"))
        processed_file <- file.path(project_dir, paste0(project_name, "_processed.rds"))
        
        # Determine status based on files
        if (file.exists(session_file) && file.exists(processed_file)) {
          status <- "Completed"
        } else if (file.exists(session_file)) {
          status <- "In Progress"
        } else {
          status <- "Unknown"
        }
        
        # Get creation date
        creation_time <- file.info(project_dir)$ctime
        date_str <- format(creation_time, "%Y-%m-%d")
        
        log_to_terminal(paste("Project:", project_name, "Status:", status, "Date:", date_str), "INFO")
        
        data.frame(
          Project = project_name,
          Date = date_str,
          Status = status,
          stringsAsFactors = FALSE
        )
      })
      
      # Combine and sort by date (most recent first)
      all_projects <- do.call(rbind, project_info)
      all_projects <- all_projects[order(all_projects$Date, decreasing = TRUE), ]
      
      # Return top 5 projects
      result <- head(all_projects, 5)
      log_to_terminal(paste("Returning", nrow(result), "projects"), "INFO")
      return(result)
    } else {
      log_to_terminal("No project directories found", "INFO")
      # No projects found
      data.frame(
        Project = character(0),
        Date = character(0),
        Status = character(0),
        stringsAsFactors = FALSE
      )
    }
  })
  
  # Pipeline execution
  observeEvent(input$run_pipeline, {
    if (values$pipeline_running) {
      showNotification("Pipeline is already running!", type = "warning")
      return()
    }
    
    # Log to terminal
    log_to_terminal("=== Starting Pipeline Execution ===", "INFO")
    
    # Validate inputs with better error handling
    validation_result <- validate_inputs()
    if (!validation_result$valid) {
      log_to_terminal(paste("Validation failed:", validation_result$message), "ERROR")
      showNotification(validation_result$message, type = "error")
      return()
    }
    
    log_to_terminal("Input validation passed", "INFO")
    
    values$pipeline_running <- TRUE
    values$status <- "running"
    values$log_messages <- character()
    values$progress <- 0
    
    # Create command arguments
    args <- create_pipeline_args()
    
    # Log arguments to terminal
    log_to_terminal(paste("Project name:", args$project_name), "INFO")
    log_to_terminal(paste("Working directory:", args$working_dir), "INFO")
    if (!is.null(args$rds_input)) {
      log_to_terminal(paste("RDS input:", args$rds_input), "INFO")
    }
    if (!is.null(args$h5_input)) {
      log_to_terminal(paste("H5 input:", args$h5_input), "INFO")
    }
    
    # Show initial notification
    showNotification("Starting pipeline...", type = "default", duration = NULL, id = "pipeline_status")
    
    # Run pipeline in background using promises
    future::plan(future::multisession)
    
    # Create a promise for the pipeline execution
    pipeline_promise <- future::future({
      # Run the actual pipeline with comprehensive error handling
      tryCatch({
        log_to_terminal("Initializing pipeline in background...", "INFO")
        
        # Add detailed logging
        isolate({
          values$log_messages <- c(values$log_messages, "Initializing pipeline...")
          values$log_messages <- c(values$log_messages, paste("Project name:", args$project_name))
          values$log_messages <- c(values$log_messages, paste("Working directory:", args$working_dir))
        })
        
        log_to_terminal("Calling multi-sample pipeline function...", "INFO")
        result <- run_multi_sample_pipeline(args)
        log_to_terminal("multi-sample pipeline function completed", "INFO")
        
        # Handle multi-sample pipeline result
        if (!is.null(result) && is.list(result)) {
          log_to_terminal("Multi-sample pipeline completed successfully", "INFO")
          
          # The result should already be in the correct format for multi-sample analysis
          if (!is.null(result$project_name) && !is.null(result$samples)) {
            log_to_terminal("Multi-sample results structure is valid", "INFO")
            result  # Return the result as-is
          } else {
            log_to_terminal("Multi-sample results structure is invalid", "ERROR")
            list(
              project_name = args$project_name,
              project_output_dir = file.path(args$working_dir, paste0(args$project_name, "_outputs")),
              samples = list(),
              summary = data.frame(),
              n_samples = 0,
              n_successful = 0,
              n_failed = 0,
              parallel_cores = 1,
              error = "Multi-sample results structure is invalid"
            )
          }
        } else {
          log_to_terminal("Pipeline returned null or invalid result", "ERROR")
          list(
            project_name = args$project_name,
            project_output_dir = file.path(args$working_dir, paste0(args$project_name, "_outputs")),
            samples = list(),
            summary = data.frame(),
            n_samples = 0,
            n_successful = 0,
            n_failed = 0,
            parallel_cores = 1,
            error = "Pipeline returned null or invalid result"
          )
        }
      }, error = function(e) {
        # Enhanced error logging to terminal
        log_to_terminal(paste("PIPELINE ERROR:", e$message), "ERROR")
        log_to_terminal(paste("Error occurred at:", Sys.time()), "ERROR")
        log_to_terminal(paste("Error call:", deparse(e$call)), "ERROR")
        log_to_terminal(paste("Working directory:", getwd()), "ERROR")
        log_to_terminal(paste("Project name:", args$project_name), "ERROR")
        
        # Also log to Shiny interface
        error_msg <- paste("Pipeline error:", e$message)
        isolate({
          values$log_messages <- c(values$log_messages, error_msg)
          values$log_messages <- c(values$log_messages, paste("Error occurred at:", Sys.time()))
          values$log_messages <- c(values$log_messages, paste("Error call:", deparse(e$call)))
        })
        
        # Return error information
        list(
          project_name = args$project_name,
          project_output_dir = file.path(args$working_dir, paste0(args$project_name, "_outputs")),
          samples = list(),
          summary = data.frame(),
          n_samples = 0,
          n_successful = 0,
          n_failed = 0,
          parallel_cores = 1,
          error = e$message
        )
      })
    })
    
    # Handle the promise result
    promises::then(pipeline_promise, 
      onFulfilled = function(result) {
        # Check if pipeline returned an error
        if (!is.null(result$error)) {
          # Handle pipeline error
          log_to_terminal(paste("Pipeline failed in promise handler:", result$error), "ERROR")
          isolate({
            values$pipeline_running <- FALSE
            values$status <- "error"
            values$progress <- 0
            values$log_messages <- c(values$log_messages, paste("Pipeline error:", result$error))
          })
          showNotification(paste("Pipeline failed:", result$error), type = "error", id = "pipeline_status")
        } else {
          # Update reactive values in the main session
          log_to_terminal("Pipeline completed successfully in promise handler", "INFO")
          isolate({
            values$results <- result
            values$pipeline_running <- FALSE
            values$status <- "completed"
            values$progress <- 100  # Complete the progress bar
            values$log_messages <- c(values$log_messages, "Pipeline completed successfully!")
          })
          
          # Show success notification
          showNotification("Pipeline completed successfully!", type = "default", id = "pipeline_status")
        }
      },
      onRejected = function(error) {
        # Handle errors
        log_to_terminal(paste("Promise rejected with error:", error$message), "ERROR")
        isolate({
          values$pipeline_running <- FALSE
          values$status <- "error"
          values$progress <- 0  # Reset progress on error
          values$log_messages <- c(values$log_messages, paste("Error:", error$message))
        })
        
        # Show error notification
        showNotification(paste("Pipeline failed:", error$message), type = "error", id = "pipeline_status")
      }
    )
  })
  
  # Pipeline log
  output$pipeline_log <- renderPrint({
    if (length(values$log_messages) > 0) {
      cat(paste(values$log_messages, collapse = "\n"))
    }
  })
  
  # Status indicator
  output$pipeline_status_indicator <- renderUI({
    status <- values$status
    if (status == "running") {
      div(
        style = "color: #007bff; font-weight: bold;",
        icon("spinner", class = "fa-spin"),
        " Pipeline is running..."
      )
    } else if (status == "completed") {
      div(
        style = "color: #28a745; font-weight: bold;",
        icon("check-circle"),
        " Pipeline completed successfully!"
      )
    } else if (status == "error") {
      div(
        style = "color: #dc3545; font-weight: bold;",
        icon("exclamation-triangle"),
        " Pipeline failed. Check logs for details."
      )
    } else {
      div(
        style = "color: #6c757d;",
        icon("info-circle"),
        " Ready to run pipeline"
      )
    }
  })
  
  # Sample preview
  output$sample_preview <- renderPrint({
    if (input$input_type == "h5" && !is.null(input$h5_files)) {
      cat("=== H5 Files Preview ===\n")
      for (i in 1:nrow(input$h5_files)) {
        file_info <- input$h5_files[i, ]
        sample_name <- tools::file_path_sans_ext(file_info$name)
        cat(paste("Sample", i, ":", sample_name, "\n"))
        cat(paste("  File:", file_info$name, "\n"))
        cat(paste("  Size:", format(file.size(file_info$datapath), units = "MB"), "\n"))
        cat("\n")
      }
    } else if (input$input_type == "rds" && !is.null(input$rds_files)) {
      cat("=== RDS Files Preview ===\n")
      for (i in 1:nrow(input$rds_files)) {
        file_info <- input$rds_files[i, ]
        sample_name <- tools::file_path_sans_ext(file_info$name)
        cat(paste("Sample", i, ":", sample_name, "\n"))
        cat(paste("  File:", file_info$name, "\n"))
        cat(paste("  Size:", format(file.size(file_info$datapath), units = "MB"), "\n"))
        cat("\n")
      }
    } else {
      cat("No files selected. Please upload files to see sample preview.\n")
    }
  })
  
  # Refresh samples button
  observeEvent(input$refresh_samples, {
    # This will trigger a re-render of the sample preview
  })
  
  # Progress bar with better tracking
  observe({
    if (values$pipeline_running) {
      # Animate progress bar
      invalidateLater(1000)  # Update every second
      current_progress <- isolate({
        if (is.null(values$progress)) values$progress <- 0
        values$progress <- min(values$progress + 10, 90)  # Increment but don't reach 100% until done
        values$progress
      })
      updateProgressBar(session, "pipeline_progress", value = current_progress)
    } else {
      updateProgressBar(session, "pipeline_progress", value = 0)
      isolate({
        values$progress <- 0
      })
    }
  })
  
  # Project summary
  output$project_summary <- renderPrint({
    if (!is.null(values$results)) {
      cat("=== Multi-Sample Pipeline Results ===\n")
      cat("Project:", values$results$project_name, "\n")
      cat("Project output directory:", values$results$project_output_dir, "\n")
      cat("Total samples:", values$results$n_samples, "\n")
      cat("Successful samples:", values$results$n_successful, "\n")
      cat("Failed samples:", values$results$n_failed, "\n")
      cat("Parallel cores used:", values$results$parallel_cores, "\n")
      cat("\n")
      
      if (!is.null(values$results$summary)) {
        cat("=== Sample Summary ===\n")
        print(values$results$summary)
      }
    } else {
      cat("No results available. Run the pipeline first.")
    }
  })
  
  # Sample selection
  observe({
    if (!is.null(values$results) && !is.null(values$results$samples)) {
      successful_samples <- names(values$results$samples)[sapply(values$results$samples, function(x) x$status == "completed")]
      updateSelectInput(session, "selected_sample", choices = successful_samples)
    }
  })
  
  # Sample summary
  output$sample_summary <- renderPrint({
    if (!is.null(input$selected_sample) && !is.null(values$results$samples)) {
      sample_result <- values$results$samples[[input$selected_sample]]
      if (!is.null(sample_result) && sample_result$status == "completed") {
        cat("=== Sample Details ===\n")
        cat("Sample name:", sample_result$sample_name, "\n")
        cat("Input file:", sample_result$input_file, "\n")
        cat("Input type:", sample_result$input_type, "\n")
        cat("Number of cells:", sample_result$n_cells, "\n")
        cat("Number of genes:", sample_result$n_genes, "\n")
        cat("Number of clusters:", sample_result$n_clusters, "\n")
        cat("Output directory:", sample_result$output_dir, "\n")
      } else {
        cat("Selected sample failed to process or not found.")
      }
    } else {
      cat("No sample selected or no results available.")
    }
  })
  
  # Load sample results
  observeEvent(input$load_sample_results, {
    if (!is.null(input$selected_sample) && !is.null(values$results$samples)) {
      sample_result <- values$results$samples[[input$selected_sample]]
      if (!is.null(sample_result) && sample_result$status == "completed") {
        # Load the sample's Seurat object for detailed analysis
        tryCatch({
          if (file.exists(sample_result$final_file)) {
            seurat_obj <- readRDS(sample_result$final_file)
            values$current_sample_seurat <- seurat_obj
            values$current_sample_name <- input$selected_sample
            showNotification(paste("Loaded sample:", input$selected_sample), type = "success")
          } else {
            showNotification(paste("Sample file not found:", sample_result$final_file), type = "error")
          }
        }, error = function(e) {
          showNotification(paste("Error loading sample:", e$message), type = "error")
        })
      }
    }
  })
  
  # Comparison plots
  output$comparison_plots <- renderUI({
    if (!is.null(values$results) && !is.null(values$results$project_output_dir)) {
      comparison_dir <- file.path(values$results$project_output_dir, "comparisons")
      if (dir.exists(comparison_dir)) {
        plot_files <- list.files(comparison_dir, pattern = "\\.png$", full.names = TRUE)
        if (length(plot_files) > 0) {
          plot_list <- lapply(plot_files, function(plot_file) {
            plot_name <- basename(plot_file)
            box(
              title = plot_name,
              width = 6,
              imageOutput(paste0("comparison_", plot_name), height = "400px")
            )
          })
          do.call(tagList, plot_list)
        } else {
          p("No comparison plots found. Run the pipeline first.")
        }
      } else {
        p("Comparison directory not found. Run the pipeline first.")
      }
    } else {
      p("No results available. Run the pipeline first.")
    }
  })
  
  # Individual comparison plot renderers
  observe({
    if (!is.null(values$results) && !is.null(values$results$project_output_dir)) {
      comparison_dir <- file.path(values$results$project_output_dir, "comparisons")
      if (dir.exists(comparison_dir)) {
        plot_files <- list.files(comparison_dir, pattern = "\\.png$", full.names = TRUE)
        lapply(plot_files, function(plot_file) {
          plot_name <- basename(plot_file)
          local({
            plot_file_local <- plot_file
            plot_name_local <- plot_name
            output[[paste0("comparison_", plot_name_local)]] <- renderImage({
              list(src = plot_file_local, width = "100%", height = "400px")
            }, deleteFile = FALSE)
          })
        })
      }
    }
  })
  
  # Sample plot outputs
  output$sample_plot_outputs <- renderUI({
    if (!is.null(values$current_sample_seurat) && !is.null(values$current_sample_name)) {
      sample_result <- values$results$samples[[values$current_sample_name]]
      if (!is.null(sample_result) && sample_result$status == "completed") {
        sample_output_dir <- sample_result$output_dir
        plot_files <- list.files(sample_output_dir, pattern = "\\.png$", full.names = TRUE)
        if (length(plot_files) > 0) {
          plot_list <- lapply(plot_files, function(plot_file) {
            plot_name <- basename(plot_file)
            box(
              title = paste(values$current_sample_name, "-", plot_name),
              width = 6,
              imageOutput(paste0("sample_", plot_name), height = "400px")
            )
          })
          do.call(tagList, plot_list)
        } else {
          p("No plots found for this sample.")
        }
      } else {
        p("Sample not processed successfully.")
      }
    } else {
      p("No sample loaded. Select and load a sample first.")
    }
  })
  
  # Individual sample plot renderers
  observe({
    if (!is.null(values$current_sample_seurat) && !is.null(values$current_sample_name)) {
      sample_result <- values$results$samples[[values$current_sample_name]]
      if (!is.null(sample_result) && sample_result$status == "completed") {
        sample_output_dir <- sample_result$output_dir
        plot_files <- list.files(sample_output_dir, pattern = "\\.png$", full.names = TRUE)
        lapply(plot_files, function(plot_file) {
          plot_name <- basename(plot_file)
          local({
            plot_file_local <- plot_file
            plot_name_local <- plot_name
            output[[paste0("sample_", plot_name_local)]] <- renderImage({
              list(src = plot_file_local, width = "100%", height = "400px")
            }, deleteFile = FALSE)
          })
        })
      }
    }
  })
  
  # Sample results table
  output$sample_results_table <- DT::renderDataTable({
    if (!is.null(values$current_sample_seurat)) {
      # Create a summary table for the current sample
      sample_summary <- data.frame(
        Metric = c("Sample Name", "Number of Cells", "Number of Genes", "Number of Clusters"),
        Value = c(
          values$current_sample_name,
          ncol(values$current_sample_seurat),
          nrow(values$current_sample_seurat),
          length(unique(Seurat::Idents(values$current_sample_seurat)))
        )
      )
      DT::datatable(sample_summary, options = list(pageLength = 10))
    } else {
      DT::datatable(data.frame(Message = "No sample loaded"))
    }
  })
  
  # Sample file outputs
  output$sample_file_outputs <- renderUI({
    if (!is.null(values$current_sample_name)) {
      sample_result <- values$results$samples[[values$current_sample_name]]
      if (!is.null(sample_result) && sample_result$status == "completed") {
        sample_output_dir <- sample_result$output_dir
        file_list <- list.files(sample_output_dir, full.names = TRUE)
        if (length(file_list) > 0) {
          file_boxes <- lapply(file_list, function(file_path) {
            file_name <- basename(file_path)
            file_size <- format(file.size(file_path), units = "MB")
            box(
              title = file_name,
              width = 6,
              p(paste("Size:", file_size)),
              downloadButton(paste0("download_", file_name), "Download")
            )
          })
          do.call(tagList, file_boxes)
        } else {
          p("No files found for this sample.")
        }
      } else {
        p("Sample not processed successfully.")
      }
    } else {
      p("No sample selected.")
    }
  })
  
  # Project file outputs
  output$project_file_outputs <- renderUI({
    if (!is.null(values$results) && !is.null(values$results$project_output_dir)) {
      project_dir <- values$results$project_output_dir
      file_list <- list.files(project_dir, full.names = TRUE, recursive = TRUE)
      if (length(file_list) > 0) {
        file_boxes <- lapply(file_list, function(file_path) {
          file_name <- basename(file_path)
          rel_path <- gsub(paste0(project_dir, "/"), "", file_path)
          file_size <- format(file.size(file_path), units = "MB")
          box(
            title = file_name,
            width = 6,
            p(paste("Path:", rel_path)),
            p(paste("Size:", file_size)),
            downloadButton(paste0("download_project_", file_name), "Download")
          )
        })
        do.call(tagList, file_boxes)
      } else {
        p("No project files found.")
      }
    } else {
      p("No results available. Run the pipeline first.")
    }
  })
  
  # Individual plot renderers
  observe({
    if (!is.null(values$results) && !is.null(values$results$plots)) {
      lapply(names(values$results$plots), function(plot_name) {
        local({
          plot_name_local <- plot_name
          output[[paste0("plot_", plot_name_local)]] <- renderImage({
            # Try to load and display the plot
            plot_file <- values$results$plots[[plot_name_local]]
            if (file.exists(plot_file)) {
              # Return the image file
              list(
                src = plot_file,
                contentType = "image/png",
                width = "100%",
                height = "auto",
                alt = plot_name_local
              )
            } else {
              # Return a placeholder image
              list(
                src = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
                contentType = "image/png",
                width = "100%",
                height = "auto",
                alt = paste("Plot not found:", plot_name_local)
              )
            }
          }, deleteFile = FALSE)
        })
      })
    }
  })
  
  # Results table
  output$results_table <- DT::renderDataTable({
    if (!is.null(values$results)) {
      # Create a summary table
      summary_data <- data.frame(
        Metric = c("Project Name", "Output Directory", "Number of Cells", "Number of Genes", "Number of Clusters"),
        Value = c(
          values$results$project_name,
          values$results$output_dir,
          as.character(values$results$n_cells),
          as.character(values$results$n_genes),
          as.character(values$results$n_clusters)
        ),
        stringsAsFactors = FALSE
      )
      DT::datatable(summary_data, options = list(pageLength = 10))
    } else {
      DT::datatable(data.frame(Message = "No results available. Run the pipeline first."))
    }
  })
  
  # File outputs with proper Shiny file serving
  output$file_outputs <- renderUI({
    if (!is.null(values$results) && !is.null(values$results$files)) {
      file_list <- lapply(values$results$files, function(file_path) {
        # Check if file exists
        if (file.exists(file_path)) {
          # Create a download link using Shiny's downloadHandler
          downloadLink(
            outputId = paste0("download_", basename(file_path)),
            label = basename(file_path)
          )
        } else {
          tags$div(
            style = "color: red;",
            paste("File not found:", basename(file_path))
          )
        }
      })
      do.call(tagList, file_list)
    } else {
      p("No files available. Run the pipeline first.")
    }
  })
  
  # Download handlers for files
  observe({
    if (!is.null(values$results) && !is.null(values$results$files)) {
      lapply(values$results$files, function(file_path) {
        if (file.exists(file_path)) {
          local({
            file_path_local <- file_path
            output[[paste0("download_", basename(file_path_local))]] <- downloadHandler(
              filename = function() {
                basename(file_path_local)
              },
              content = function(file) {
                file.copy(file_path_local, file)
              }
            )
          })
        }
      })
    }
  })
  
  # Helper functions
  validate_inputs <- function() {
    log_to_terminal("=== Validating Multi-Sample Inputs ===", "INFO")
    
    # Enhanced validation with detailed error messages for multiple files
    if (input$input_type == "h5") {
      log_to_terminal("Validating H5 inputs", "INFO")
      if (is.null(input$h5_files) || nrow(input$h5_files) == 0) {
        log_to_terminal("H5 files are null or empty", "ERROR")
        return(list(valid = FALSE, message = "Please select at least one H5 file!"))
      }
      
      log_to_terminal(paste("Number of H5 files:", nrow(input$h5_files)), "INFO")
      
      # Validate each file
      for (i in 1:nrow(input$h5_files)) {
        file_info <- input$h5_files[i, ]
        log_to_terminal(paste("Validating H5 file", i, ":", file_info$name), "INFO")
        
        # Check if file exists and is readable
        if (!file.exists(file_info$datapath)) {
          log_to_terminal(paste("H5 file not found:", file_info$datapath), "ERROR")
          return(list(valid = FALSE, message = paste("H5 file not found:", file_info$name)))
        }
        
        # Check file size
        file_size <- file.size(file_info$datapath)
        log_to_terminal(paste("H5 file size:", file_size), "INFO")
        if (file_size == 0) {
          log_to_terminal("H5 file is empty", "ERROR")
          return(list(valid = FALSE, message = paste("Selected H5 file is empty:", file_info$name)))
        }
      }
      
      log_to_terminal("All H5 files validation passed", "INFO")
    }
    
    if (input$input_type == "rds") {
      log_to_terminal("Validating RDS inputs", "INFO")
      if (is.null(input$rds_files) || nrow(input$rds_files) == 0) {
        log_to_terminal("RDS files are null or empty", "ERROR")
        return(list(valid = FALSE, message = "Please select at least one RDS file!"))
      }
      
      log_to_terminal(paste("Number of RDS files:", nrow(input$rds_files)), "INFO")
      
      # Validate each file
      for (i in 1:nrow(input$rds_files)) {
        file_info <- input$rds_files[i, ]
        log_to_terminal(paste("Validating RDS file", i, ":", file_info$name), "INFO")
        
        # Check if file exists and is readable
        if (!file.exists(file_info$datapath)) {
          log_to_terminal(paste("RDS file not found:", file_info$datapath), "ERROR")
          return(list(valid = FALSE, message = paste("RDS file not found:", file_info$name)))
        }
        
        # Check file size
        file_size <- file.size(file_info$datapath)
        log_to_terminal(paste("RDS file size:", file_size), "INFO")
        if (file_size == 0) {
          log_to_terminal("RDS file is empty", "ERROR")
          return(list(valid = FALSE, message = paste("Selected RDS file is empty:", file_info$name)))
        }
      }
      
      log_to_terminal("All RDS files validation passed", "INFO")
    }
    
    # Validate SouporCell files if provided
    if (input$use_souporcell) {
      if (is.null(input$souporcell_files) || nrow(input$souporcell_files) == 0) {
        log_to_terminal("SouporCell files are required but not provided", "ERROR")
        return(list(valid = FALSE, message = "Please select SouporCell files!"))
      }
      
      # Check if number of SouporCell files matches input files
      n_input_files <- if (input$input_type == "h5") nrow(input$h5_files) else nrow(input$rds_files)
      n_souporcell_files <- nrow(input$souporcell_files)
      
      if (n_souporcell_files != n_input_files) {
        log_to_terminal(paste("Mismatch: input files =", n_input_files, ", SouporCell files =", n_souporcell_files), "ERROR")
        return(list(valid = FALSE, message = paste("Number of SouporCell files (", n_souporcell_files, ") must match number of input files (", n_input_files, ")!")))
      }
    }
    
    # Validate sample configuration file if provided
    if (input$use_sample_configs) {
      if (is.null(input$sample_configs_file)) {
        log_to_terminal("Sample configuration file is required but not provided", "ERROR")
        return(list(valid = FALSE, message = "Please select a sample configuration file!"))
      }
      
      if (!file.exists(input$sample_configs_file$datapath)) {
        log_to_terminal(paste("Sample configuration file not found:", input$sample_configs_file$datapath), "ERROR")
        return(list(valid = FALSE, message = paste("Sample configuration file not found:", input$sample_configs_file$name)))
      }
    }
    
    # Validate project name
    log_to_terminal(paste("Project name:", input$project_name), "INFO")
    if (is.null(input$project_name) || input$project_name == "") {
      log_to_terminal("Project name is empty", "ERROR")
      return(list(valid = FALSE, message = "Please provide a project name!"))
    }
    
    # Check for invalid characters in project name
    if (grepl("[^a-zA-Z0-9_-]", input$project_name)) {
      log_to_terminal("Project name contains invalid characters", "ERROR")
      return(list(valid = FALSE, message = "Project name can only contain letters, numbers, underscores, and hyphens!"))
    }
    
    log_to_terminal("All multi-sample input validation passed", "INFO")
    return(list(valid = TRUE, message = "Validation passed"))
  }
  
  create_pipeline_args <- function() {
    log_to_terminal("=== Creating Multi-Sample Pipeline Arguments ===", "INFO")
    
    args <- list()
    
    # Input files with proper path handling for multiple files
    if (input$input_type == "h5") {
      args$h5_inputs <- sapply(1:nrow(input$h5_files), function(i) {
        normalizePath(input$h5_files$datapath[i], mustWork = FALSE)
      })
      log_to_terminal(paste("H5 input paths:", paste(args$h5_inputs, collapse = ", ")), "INFO")
    } else {
      args$rds_inputs <- sapply(1:nrow(input$rds_files), function(i) {
        normalizePath(input$rds_files$datapath[i], mustWork = FALSE)
      })
      log_to_terminal(paste("RDS input paths:", paste(args$rds_inputs, collapse = ", ")), "INFO")
    }
    
    # SouporCell files
    if (input$use_souporcell) {
      args$soupor_cell_doublet_inputs <- sapply(1:nrow(input$souporcell_files), function(i) {
        normalizePath(input$souporcell_files$datapath[i], mustWork = FALSE)
      })
      log_to_terminal(paste("SouporCell input paths:", paste(args$soupor_cell_doublet_inputs, collapse = ", ")), "INFO")
    }
    
    # Sample configuration file
    if (input$use_sample_configs) {
      args$sample_configs <- normalizePath(input$sample_configs_file$datapath, mustWork = FALSE)
      log_to_terminal(paste("Sample config path:", args$sample_configs), "INFO")
    }
    
    # Parallel processing settings
    args$parallel <- input$use_parallel
    if (input$use_parallel) {
      args$n_cores <- input$n_cores
      log_to_terminal(paste("Parallel processing enabled with", args$n_cores, "cores"), "INFO")
    }
    
    # Project settings
    args$project_name <- input$project_name
    
    # Ensure working directory is a proper absolute path
    if (is.null(input$working_dir) || input$working_dir == "" || input$working_dir == ".") {
      args$working_dir <- getwd()
    } else {
      args$working_dir <- normalizePath(input$working_dir, mustWork = FALSE)
    }
    
    # Ensure the working directory exists
    if (!dir.exists(args$working_dir)) {
      dir.create(args$working_dir, recursive = TRUE)
      log_to_terminal(paste("Created working directory:", args$working_dir), "INFO")
    }
    
    log_to_terminal(paste("Project name:", args$project_name), "INFO")
    log_to_terminal(paste("Working directory:", args$working_dir), "INFO")
    
    # QC parameters
    args$min_features <- input$min_features
    args$min_counts <- input$min_counts
    args$max_features <- input$max_features
    args$max_counts <- input$max_counts
    args$max_mt_percent <- input$max_mt_percent
    args$min_cells <- input$min_cells
    log_to_terminal(paste("Min features:", args$min_features), "INFO")
    log_to_terminal(paste("Min counts:", args$min_counts), "INFO")
    log_to_terminal(paste("Max features:", args$max_features), "INFO")
    log_to_terminal(paste("Max counts:", args$max_counts), "INFO")
    log_to_terminal(paste("Max MT percent:", args$max_mt_percent), "INFO")
    log_to_terminal(paste("Min cells:", args$min_cells), "INFO")
    
    # Processing parameters
    args$n_variable_features <- input$n_variable_features
    args$normalization_method <- input$normalization_method
    args$scaling_method <- input$scaling_method
    args$pca_dimensions <- input$pca_dimensions
    args$scale_factor <- input$scale_factor
    log_to_terminal(paste("Variable features:", args$n_variable_features), "INFO")
    log_to_terminal(paste("Normalization method:", args$normalization_method), "INFO")
    log_to_terminal(paste("Scaling method:", args$scaling_method), "INFO")
    log_to_terminal(paste("PCA dimensions:", args$pca_dimensions), "INFO")
    log_to_terminal(paste("Scale factor:", args$scale_factor), "INFO")
    
    # Clustering parameters
    args$clustering_resolution <- input$clustering_resolution
    args$clustering_algorithm <- input$clustering_algorithm
    args$min_cluster_size <- input$min_cluster_size
    log_to_terminal(paste("Clustering resolution:", args$clustering_resolution), "INFO")
    log_to_terminal(paste("Clustering algorithm:", args$clustering_algorithm), "INFO")
    log_to_terminal(paste("Min cluster size:", args$min_cluster_size), "INFO")
    
    # UMAP parameters
    args$umap_n_neighbors <- input$umap_n_neighbors
    args$umap_min_dist <- input$umap_min_dist
    log_to_terminal(paste("UMAP neighbors:", args$umap_n_neighbors), "INFO")
    log_to_terminal(paste("UMAP min dist:", args$umap_min_dist), "INFO")
    
    # Analysis options
    args$find_markers <- input$find_markers
    args$save_intermediate <- input$save_intermediate
    args$verbose <- input$verbose
    log_to_terminal(paste("Find markers:", args$find_markers), "INFO")
    log_to_terminal(paste("Save intermediate:", args$save_intermediate), "INFO")
    log_to_terminal(paste("Verbose:", args$verbose), "INFO")
    
    # Configuration file
    config_file_path <- "config/settings.yaml"
    if (file.exists(config_file_path)) {
      args$config_file <- normalizePath(config_file_path, mustWork = FALSE)
      log_to_terminal(paste("Config file:", args$config_file), "INFO")
    } else {
      log_to_terminal("Warning: Config file not found, using defaults", "WARN")
      args$config_file <- NULL
    }
    
    log_to_terminal("Pipeline arguments created successfully", "INFO")
    return(args)
  }
  

}

# Run the application
shinyApp(ui = ui, server = server)

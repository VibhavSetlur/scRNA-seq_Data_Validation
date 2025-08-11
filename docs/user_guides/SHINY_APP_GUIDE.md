# Shiny App User Guide

## What You'll See When You Open the App

The snRNA-seq Pipeline Shiny app provides a comprehensive, user-friendly interface for single-nucleus RNA sequencing analysis. Here's what you'll find:

## 🏠 Dashboard Tab

**Welcome Message**: Introduction to the pipeline and its capabilities

**System Status Panel**: 
- R version and platform information
- Available memory and CPU cores
- Working directory status

**Recent Projects Panel**: 
- List of previously run analyses
- Project status (Completed, Running, Failed)

## 📁 Data Input Tab

**Input Type Selection**:
- **H5 File**: For 10X Genomics data (.h5, .hdf5 files)
- **RDS File**: For pre-processed Seurat objects (.rds files)

**File Upload**:
- Drag-and-drop or click-to-browse file selection
- File validation and preview

**Project Settings**:
- Project name input
- Working directory selection
- SouporCell doublet detection (optional)

## 🔍 Quality Control Tab

**QC Parameters**:
- **Minimum Features per Cell**: 200 (default)
- **Minimum Counts per Cell**: 1000 (default)
- **Maximum Features per Cell**: 6000 (default)
- **Maximum Counts per Cell**: 25000 (default)
- **Maximum Mitochondrial %**: 20% (default)
- **Minimum Cells per Gene**: 3 (default)

**Real-time Validation**: Parameters are validated as you type

## ⚙️ Processing Tab

**Normalization Options**:
- **CLR**: Centered Log Ratio (recommended for UMI data)
- **LogNormalize**: Log normalization
- **RC**: Relative Counts

**Feature Selection**:
- **Number of Variable Features**: 2000 (default)
- **Scaling Method**: negbinom (recommended) or linear
- **PCA Dimensions**: 15 (default)
- **Scale Factor**: 10000 (default)

**Output Options**:
- Save intermediate results toggle

## 🎯 Clustering Tab

**Clustering Algorithm**:
- **leiden**: Recommended for most datasets
- **louvain**: Alternative option
- **multilevel**: For specific cases
- **slm**: For large datasets

**Clustering Parameters**:
- **Resolution**: 0.5 (default, controls number of clusters)
- **Minimum Cluster Size**: 10 (default)
- **UMAP n_neighbors**: 30 (default)
- **UMAP min_dist**: 0.3 (default)
- **Random Seed**: 42 (default)

## 📊 Visualization Tab

**Plot Settings**:
- **Theme**: classic, minimal, or bw
- **Color Palette**: viridis (recommended), rainbow, or brewer
- **Point Size**: 0.7 (default)
- **Label Size**: 3 (default)
- **Title Size**: 14 (default)

**Output Format**:
- **PNG**: Portable Network Graphics
- **PDF**: Portable Document Format
- **SVG**: Scalable Vector Graphics

## ▶️ Run Pipeline Tab

**Analysis Options**:
- **Find Cluster Markers**: Toggle for differential expression analysis
- **Verbose Output**: Toggle for detailed progress information

**Execution**:
- Large "Run Pipeline" button
- Real-time progress bar
- Live log output
- Status indicators

## 📈 Results Tab

**Summary Panel**:
- Project completion status
- Number of cells and genes
- Number of clusters identified
- Output directory location

**Plots Panel**:
- Interactive plot viewer
- Download options for all plots
- Plot descriptions and interpretations

**Tables Panel**:
- Interactive data tables
- Sortable and searchable results
- Export functionality

**Files Panel**:
- Download links for all output files
- File descriptions and formats
- Direct access to results

## ❓ Help Tab

**Getting Started**: Step-by-step guide for new users

**Parameters**: Detailed explanations of all parameters

**Troubleshooting**: Common issues and solutions

**FAQ**: Frequently asked questions

## 🎨 User Interface Features

**Responsive Design**: Works on desktop, tablet, and mobile

**Real-time Validation**: Parameters are checked as you enter them

**Progress Tracking**: Live updates during pipeline execution

**Error Handling**: User-friendly error messages and suggestions

**Download Options**: Easy access to all results and plots

## 🚀 Getting Started

1. **Upload Your Data**: Choose H5 or RDS file format
2. **Configure Parameters**: Use defaults or adjust based on your data
3. **Run Analysis**: Click "Run Pipeline" and monitor progress
4. **View Results**: Explore plots, tables, and download files

## 💡 Tips for Best Results

- **Start with Defaults**: Default parameters work well for most datasets
- **Check QC Plots**: Use quality control visualizations to guide parameter selection
- **Monitor Progress**: Watch the pipeline log for any issues
- **Save Results**: Download and save your results for future reference

## 🔧 Troubleshooting

If you encounter issues:
1. Check the Help tab for guidance
2. Review the troubleshooting section
3. Check the pipeline log for error messages
4. Try with more conservative parameters

The Shiny app is designed to be intuitive and user-friendly while providing access to all the powerful features of the snRNA-seq pipeline!

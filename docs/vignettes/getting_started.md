# Getting Started Guide

## Quick Start

This guide will help you get up and running with the snRNA-seq Pipeline in minutes.

### 1. Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/snRNA-seq-Pipeline.git
cd snRNA-seq-Pipeline

# Run the setup script
chmod +x setup/setup.sh
./setup/setup.sh
```

### 2. Test Your Installation

```bash
# Run the test suite
./tests/run_tests.sh
```

### 3. Run Your First Analysis

```bash
# Basic analysis with test data
Rscript src/core/pipeline.R \
  --h5_input tests/data/test_data.h5 \
  --project_name MyFirstAnalysis \
  --working_dir results/
```

## Next Steps

- Read the [Parameter Tuning Guide](parameter_tuning.md)
- Explore the [Shiny Application](shiny_app.md)
- Check out [Advanced Analysis](advanced_analysis.md)

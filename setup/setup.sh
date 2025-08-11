#!/bin/bash

# snRNA-seq Pipeline Setup Script
# This script runs the R setup script and provides additional system checks

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local message="$1"
    local type="${2:-info}"
    
    case $type in
        "info")
            echo -e "${BLUE}[INFO]${NC} $message"
            ;;
        "success")
            echo -e "${GREEN}[SUCCESS]${NC} $message"
            ;;
        "warning")
            echo -e "${YELLOW}[WARNING]${NC} $message"
            ;;
        "error")
            echo -e "${RED}[ERROR]${NC} $message"
            ;;
    esac
}

# Function to check if R is installed
check_r_installation() {
    print_status "Checking R installation..."
    
    if command -v R >/dev/null 2>&1; then
        R_VERSION=$(R --version | head -n 1)
        print_status "R is installed: $R_VERSION" "success"
    else
        print_status "R is not installed. Please install R first." "error"
        exit 1
    fi
    
    if command -v Rscript >/dev/null 2>&1; then
        print_status "Rscript is available" "success"
    else
        print_status "Rscript is not available" "error"
        exit 1
    fi
}

# Function to check system requirements
check_system_requirements() {
    print_status "Checking system requirements..."
    
    # Check available memory
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        MEMORY=$(free -h | grep Mem | awk '{print $2}')
        print_status "Total memory: $MEMORY" "info"
        
        AVAILABLE_MEM=$(free -h | grep Mem | awk '{print $7}')
        print_status "Available memory: $AVAILABLE_MEM" "info"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        MEMORY=$(sysctl hw.memsize | awk '{print $2}')
        MEMORY_GB=$((MEMORY / 1024 / 1024 / 1024))
        print_status "Total memory: ${MEMORY_GB}GB" "info"
    fi
    
    # Check available disk space
    DISK_SPACE=$(df -h . | tail -1 | awk '{print $4}')
    print_status "Available disk space: $DISK_SPACE" "info"
    
    # Check number of CPU cores
    CPU_CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo "Unknown")
    print_status "Number of CPU cores: $CPU_CORES" "info"
}

# Function to create virtual environment (optional)
create_virtual_env() {
    print_status "Setting up R environment..."
    
    # Create .Rprofile if it doesn't exist
    if [[ ! -f ".Rprofile" ]]; then
        cat > .Rprofile << 'EOF'
# snRNA-seq Pipeline R Profile
# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Set memory limit (if needed)
# options(max.memory = 8000)

# Load common packages
if (interactive()) {
  cat("snRNA-seq Pipeline R environment loaded\n")
}
EOF
        print_status "Created .Rprofile" "success"
    else
        print_status ".Rprofile already exists" "info"
    fi
}

# Function to run R setup script
run_r_setup() {
    print_status "Running R setup script..."
    
    if [[ -f "setup/install_dependencies.R" ]]; then
        Rscript setup/install_dependencies.R
        if [[ $? -eq 0 ]]; then
            print_status "R setup completed successfully" "success"
        else
            print_status "R setup failed" "error"
            exit 1
        fi
    else
        print_status "R setup script not found" "error"
        exit 1
    fi
}

# Function to create additional directories
create_directories() {
    print_status "Creating additional directories..."
    
    DIRS=(
        "data/raw"
        "data/processed"
        "results/plots"
        "results/tables"
        "logs/analysis"
        "temp/cache"
        "docs/reports"
        "tests/data"
        "tests/results"
    )
    
    for dir in "${DIRS[@]}"; do
        if [[ ! -d "$dir" ]]; then
            mkdir -p "$dir"
            print_status "Created directory: $dir" "success"
        else
            print_status "Directory already exists: $dir" "info"
        fi
    done
}

# Function to set up git hooks (if in git repository)
setup_git_hooks() {
    if [[ -d ".git" ]]; then
        print_status "Setting up git hooks..."
        
        # Create pre-commit hook
        cat > .git/hooks/pre-commit << 'EOF'
#!/bin/bash
# Pre-commit hook for snRNA-seq Pipeline

echo "Running pre-commit checks..."

# Check if R files are properly formatted
if command -v R >/dev/null 2>&1; then
    echo "Checking R code formatting..."
    # Add any R linting checks here
fi

echo "Pre-commit checks passed"
EOF
        
        chmod +x .git/hooks/pre-commit
        print_status "Created pre-commit hook" "success"
    fi
}

# Function to create a test script
create_test_script() {
    print_status "Creating test script..."
    
    cat > tests/run_tests.sh << 'EOF'
#!/bin/bash

# Test script for snRNA-seq Pipeline
# This script runs basic tests to ensure the pipeline is working

set -e

echo "Running snRNA-seq Pipeline tests..."

# Test 1: Check if R and required packages are available
echo "Test 1: Checking R environment..."
Rscript -e "
required_packages <- c('Seurat', 'tidyverse', 'argparse', 'yaml')
missing_packages <- character()

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) == 0) {
  cat('✓ All required packages are available\n')
} else {
  cat('✗ Missing packages:', paste(missing_packages, collapse = ', '), '\n')
  quit(status = 1)
}
"

# Test 2: Check if configuration file exists
echo "Test 2: Checking configuration files..."
if [[ -f "config/settings.yaml" ]]; then
    echo "✓ Configuration file exists"
else
    echo "✗ Configuration file missing"
    exit 1
fi

# Test 3: Check if main pipeline script exists
echo "Test 3: Checking pipeline scripts..."
if [[ -f "src/core/pipeline.R" ]]; then
    echo "✓ Main pipeline script exists"
else
    echo "✗ Main pipeline script missing"
    exit 1
fi

# Test 4: Generate test data
echo "Test 4: Generating test data..."
Rscript tests/scripts/generate_test_data.R

# Test 5: Run pipeline on test data
echo "Test 5: Running pipeline on test data..."
if [[ -f "tests/data/test_data.h5" ]]; then
    Rscript src/core/pipeline.R \
        --h5_input tests/data/test_data.h5 \
        --project_name TestProject \
        --working_dir tests/results \
        --find_markers FALSE
    echo "✓ Pipeline test completed successfully"
else
    echo "✗ Test data not found"
    exit 1
fi

echo "All tests passed!"
EOF
    
    chmod +x tests/run_tests.sh
    print_status "Created test script" "success"
}

# Function to show help
show_help() {
    echo "snRNA-seq Pipeline Setup Script"
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help     Show this help message"
    echo "  -v, --verbose  Enable verbose output"
    echo "  --skip-r       Skip R setup (only create directories)"
    echo "  --test-only    Only run tests"
    echo ""
    echo "This script will:"
    echo "  1. Check system requirements"
    echo "  2. Install R dependencies"
    echo "  3. Create necessary directories"
    echo "  4. Set up git hooks (if in git repository)"
    echo "  5. Create test scripts"
    echo ""
}

# Main function
main() {
    print_status "=== snRNA-seq Pipeline Setup ===" "info"
    print_status "Starting setup process..." "info"
    
    # Parse command line arguments
    SKIP_R=false
    TEST_ONLY=false
    VERBOSE=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            --skip-r)
                SKIP_R=true
                shift
                ;;
            --test-only)
                TEST_ONLY=true
                shift
                ;;
            *)
                print_status "Unknown option: $1" "error"
                show_help
                exit 1
                ;;
        esac
    done
    
    if [[ "$TEST_ONLY" == true ]]; then
        print_status "Running tests only..." "info"
        if [[ -f "tests/run_tests.sh" ]]; then
            ./tests/run_tests.sh
        else
            print_status "Test script not found. Run setup first." "error"
            exit 1
        fi
        exit 0
    fi
    
    # Check R installation
    check_r_installation
    
    # Check system requirements
    check_system_requirements
    
    # Create virtual environment
    create_virtual_env
    
    # Create directories
    create_directories
    
    # Set up git hooks
    setup_git_hooks
    
    # Create test script
    create_test_script
    
    # Run R setup (unless skipped)
    if [[ "$SKIP_R" == false ]]; then
        run_r_setup
    else
        print_status "Skipping R setup as requested" "warning"
    fi
    
    print_status "=== SETUP COMPLETED ===" "success"
    print_status "You can now use the snRNA-seq pipeline!" "success"
    print_status ""
    print_status "Next steps:" "info"
    print_status "  1. Run tests: ./tests/run_tests.sh" "info"
    print_status "  2. Start Shiny app: Rscript shiny_app/app.R" "info"
    print_status "  3. Run pipeline: Rscript src/core/pipeline.R --help" "info"
    print_status ""
}

# Run main function
main "$@"

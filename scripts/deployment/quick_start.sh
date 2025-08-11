#!/bin/bash

# Quick Start Script for snRNA-seq Pipeline
# This script provides a seamless setup experience for users

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== snRNA-seq Pipeline Quick Start ===${NC}"
echo ""

# Check if R is installed
if ! command -v R &> /dev/null; then
    echo -e "${RED}Error: R is not installed. Please install R first.${NC}"
    echo "Visit: https://cran.r-project.org/"
    exit 1
fi

echo -e "${GREEN}✓ R is installed${NC}"

# Check if Docker is available (optional)
if command -v docker &> /dev/null; then
    echo -e "${GREEN}✓ Docker is available${NC}"
    DOCKER_AVAILABLE=true
else
    echo -e "${YELLOW}⚠ Docker not found (optional for containerized deployment)${NC}"
    DOCKER_AVAILABLE=false
fi

# Create necessary directories
echo -e "${BLUE}Creating directories...${NC}"
mkdir -p data results logs temp
echo -e "${GREEN}✓ Directories created${NC}"

# Install R packages if needed
echo -e "${BLUE}Checking R packages...${NC}"
if Rscript -e "library(argparse)" 2>/dev/null; then
    echo -e "${GREEN}✓ Required packages are installed${NC}"
else
    echo -e "${YELLOW}Installing required packages...${NC}"
    Rscript -e "install.packages('argparse', lib='~/.local/lib/R/library')"
    echo -e "${GREEN}✓ Packages installed${NC}"
fi

echo ""
echo -e "${GREEN}=== Setup Complete! ===${NC}"
echo ""
echo -e "${BLUE}Available options:${NC}"
echo "1. ${GREEN}Run Shiny App${NC}: ./run_shiny_app.R"
echo "2. ${GREEN}Run Terminal Pipeline${NC}: ./run_pipeline_terminal.R --help"
echo "3. ${GREEN}Docker Deployment${NC}: ./deploy.sh deploy (if Docker available)"
echo ""
echo -e "${BLUE}Quick Commands:${NC}"
echo "• Start Shiny app: Rscript run_shiny_app.R"
echo "• View help: Rscript run_pipeline_terminal.R --help"
echo "• Deploy with Docker: ./deploy.sh deploy"
echo ""
echo -e "${YELLOW}The Shiny app will be available at: http://localhost:3838${NC}"
echo ""

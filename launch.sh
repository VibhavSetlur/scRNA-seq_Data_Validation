#!/bin/bash

# Main Launcher for snRNA-seq Pipeline
# This script provides easy access to all pipeline functionality

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to show usage
show_usage() {
    echo -e "${BLUE}=== snRNA-seq Pipeline Launcher ===${NC}"
    echo ""
    echo -e "${GREEN}Available Commands:${NC}"
    echo ""
    echo -e "${BLUE}Setup & Installation:${NC}"
    echo "  ${GREEN}setup${NC}     - Run quick setup and install dependencies"
    echo "  ${GREEN}deploy${NC}    - Deploy with Docker (production ready)"
    echo ""
    echo -e "${BLUE}Running the Pipeline:${NC}"
    echo "  ${GREEN}shiny${NC}     - Launch Shiny web interface"
    echo "  ${GREEN}terminal${NC}  - Run pipeline from command line"
    echo "  ${GREEN}help${NC}      - Show command line help"
    echo ""
    echo -e "${BLUE}Examples:${NC}"
    echo "  ./launch.sh setup     # First time setup"
    echo "  ./launch.sh shiny     # Start web interface"
    echo "  ./launch.sh deploy    # Deploy with Docker"
    echo ""
}

# Function to run setup
run_setup() {
    echo -e "${BLUE}Running setup...${NC}"
    ./scripts/deployment/quick_start.sh
}

# Function to run deployment
run_deploy() {
    echo -e "${BLUE}Running deployment...${NC}"
    ./scripts/deployment/deploy.sh deploy
}

# Function to launch Shiny app
run_shiny() {
    echo -e "${BLUE}Launching Shiny app...${NC}"
    ./scripts/deployment/launch_app.sh
}

# Function to run terminal pipeline
run_terminal() {
    echo -e "${BLUE}Running terminal pipeline...${NC}"
    Rscript scripts/run_pipeline_terminal.R --help
}

# Function to show help
run_help() {
    Rscript scripts/run_pipeline_terminal.R --help
}

# Main script logic
case "${1:-}" in
    setup)
        run_setup
        ;;
    deploy)
        run_deploy
        ;;
    shiny)
        run_shiny
        ;;
    terminal)
        run_terminal
        ;;
    help)
        run_help
        ;;
    *)
        show_usage
        ;;
esac

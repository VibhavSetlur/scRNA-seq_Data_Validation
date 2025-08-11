#!/bin/bash

# snRNA-seq Pipeline Deployment Script
# This script automates the deployment of the pipeline on servers

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
DEFAULT_PORT=3838
DEFAULT_HOST="0.0.0.0"
DOCKER_IMAGE_NAME="snrna-pipeline"
CONTAINER_NAME="snrna-pipeline"

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check system requirements
check_requirements() {
    print_status "Checking system requirements..."
    
    # Check Docker
    if ! command_exists docker; then
        print_error "Docker is not installed. Please install Docker first."
        exit 1
    fi
    
    # Check Docker Compose
    if ! command_exists docker-compose; then
        print_error "Docker Compose is not installed. Please install Docker Compose first."
        exit 1
    fi
    
    # Check if Docker daemon is running
    if ! docker info >/dev/null 2>&1; then
        print_error "Docker daemon is not running. Please start Docker first."
        exit 1
    fi
    
    print_success "System requirements met"
}

# Function to create necessary directories
create_directories() {
    print_status "Creating necessary directories..."
    
    mkdir -p data
    mkdir -p results
    mkdir -p logs
    mkdir -p temp
    mkdir -p config
    mkdir -p ssl
    
    print_success "Directories created"
}

# Function to build Docker image
build_image() {
    print_status "Building Docker image..."
    
    docker build -t $DOCKER_IMAGE_NAME .
    
    if [ $? -eq 0 ]; then
        print_success "Docker image built successfully"
    else
        print_error "Failed to build Docker image"
        exit 1
    fi
}

# Function to start services
start_services() {
    local mode=$1
    
    print_status "Starting services in $mode mode..."
    
    if [ "$mode" = "production" ]; then
        docker-compose --profile production up -d
    else
        docker-compose up -d
    fi
    
    if [ $? -eq 0 ]; then
        print_success "Services started successfully"
    else
        print_error "Failed to start services"
        exit 1
    fi
}

# Function to check service health
check_health() {
    print_status "Checking service health..."
    
    # Wait for services to start
    sleep 10
    
    # Check if container is running
    if docker ps | grep -q $CONTAINER_NAME; then
        print_success "Container is running"
    else
        print_error "Container is not running"
        exit 1
    fi
    
    # Check if service is responding
    if curl -f http://localhost:$DEFAULT_PORT >/dev/null 2>&1; then
        print_success "Service is responding"
    else
        print_warning "Service is not responding yet. It may take a few minutes to start."
    fi
}

# Function to show usage information
show_usage() {
    echo "Usage: $0 [OPTIONS] COMMAND"
    echo ""
    echo "Commands:"
    echo "  build     Build the Docker image"
    echo "  start     Start the services"
    echo "  stop      Stop the services"
    echo "  restart   Restart the services"
    echo "  logs      Show service logs"
    echo "  status    Show service status"
    echo "  deploy    Full deployment (build + start)"
    echo "  clean     Remove containers and images"
    echo ""
    echo "Options:"
    echo "  -m, --mode MODE    Deployment mode (development|production) [default: development]"
    echo "  -p, --port PORT    Port number [default: $DEFAULT_PORT]"
    echo "  -h, --help         Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 deploy                    # Deploy in development mode"
    echo "  $0 deploy -m production      # Deploy in production mode"
    echo "  $0 start -p 8080             # Start on port 8080"
}

# Function to stop services
stop_services() {
    print_status "Stopping services..."
    docker-compose down
    print_success "Services stopped"
}

# Function to show logs
show_logs() {
    print_status "Showing service logs..."
    docker-compose logs -f
}

# Function to show status
show_status() {
    print_status "Service status:"
    docker-compose ps
    
    echo ""
    print_status "Container logs (last 10 lines):"
    docker-compose logs --tail=10
}

# Function to clean up
clean_up() {
    print_status "Cleaning up containers and images..."
    docker-compose down --rmi all --volumes --remove-orphans
    print_success "Cleanup completed"
}

# Main script logic
main() {
    local command=""
    local mode="development"
    local port=$DEFAULT_PORT
    
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -m|--mode)
                mode="$2"
                shift 2
                ;;
            -p|--port)
                port="$2"
                shift 2
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            build|start|stop|restart|logs|status|deploy|clean)
                command="$1"
                shift
                ;;
            *)
                print_error "Unknown option: $1"
                show_usage
                exit 1
                ;;
        esac
    done
    
    # Set environment variables
    export SHINY_PORT=$port
    export SHINY_HOST=$DEFAULT_HOST
    
    # Execute command
    case $command in
        build)
            check_requirements
            build_image
            ;;
        start)
            check_requirements
            create_directories
            start_services $mode
            check_health
            ;;
        stop)
            stop_services
            ;;
        restart)
            stop_services
            start_services $mode
            check_health
            ;;
        logs)
            show_logs
            ;;
        status)
            show_status
            ;;
        deploy)
            check_requirements
            create_directories
            build_image
            start_services $mode
            check_health
            print_success "Deployment completed!"
            echo ""
            echo "Access the application at: http://localhost:$port"
            echo "For production deployment, configure your domain and SSL certificates."
            ;;
        clean)
            clean_up
            ;;
        "")
            print_error "No command specified"
            show_usage
            exit 1
            ;;
        *)
            print_error "Unknown command: $command"
            show_usage
            exit 1
            ;;
    esac
}

# Run main function
main "$@"

# Deployment Guide - snRNA-seq Pipeline

## Overview

This guide provides comprehensive instructions for deploying the snRNA-seq Pipeline in various environments, from local development to production servers. The pipeline has been scaled and optimized to work without localhost dependencies and can be deployed on any server infrastructure.

## What Was Fixed

### 1. Shiny App Issues
- **Problem**: `object 'server' not found` error
- **Solution**: Fixed source paths in `run_shiny_app.R` to properly load the app.R file
- **Problem**: Missing markdown documentation files
- **Solution**: Created comprehensive documentation files in `shiny_app/docs/`

### 2. Server Deployment Issues
- **Problem**: Localhost dependency preventing server deployment
- **Solution**: Configured Shiny app to use `0.0.0.0` host and environment variables
- **Problem**: Missing containerization
- **Solution**: Created Docker configuration with proper networking

### 3. Terminal Interface Issues
- **Problem**: Missing argparse package
- **Solution**: Added proper library path configuration and package installation
- **Problem**: No standalone terminal interface
- **Solution**: Created `run_pipeline_terminal.R` with full command-line interface

## Deployment Options

### Option 1: Docker Deployment (Recommended)

#### Quick Start
```bash
# Clone the repository
git clone <repository-url>
cd snRNA-seq-Pipeline

# Deploy using automated script
./deploy.sh deploy
```

#### Advanced Deployment
```bash
# Development mode
./deploy.sh deploy -m development

# Production mode with NGINX
./deploy.sh deploy -m production

# Custom port
./deploy.sh deploy -p 8080
```

#### Docker Commands
```bash
# Build image
docker build -t snrna-pipeline .

# Run container
docker run -p 3838:3838 -v $(pwd)/data:/app/data snrna-pipeline

# Use docker-compose
docker-compose up -d
```

### Option 2: Local Installation

#### Prerequisites
- R 4.3.0 or higher
- Required system libraries (see Dockerfile for list)

#### Installation
```bash
# Install dependencies
Rscript setup/install_dependencies.R

# Run Shiny app
Rscript run_shiny_app.R

# Run terminal interface
Rscript run_pipeline_terminal.R --help
```

### Option 3: Server Deployment

#### Manual Server Setup
```bash
# Install system dependencies
sudo apt update
sudo apt install r-base r-base-dev gdebi-core

# Install Shiny Server
wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb
sudo gdebi shiny-server-1.5.20.1002-amd64.deb

# Deploy application
sudo cp -r shiny_app /srv/shiny-server/snrna-pipeline
sudo systemctl restart shiny-server
```

#### NGINX Reverse Proxy
```bash
# Install NGINX
sudo apt install nginx

# Configure reverse proxy
sudo cp nginx.conf /etc/nginx/sites-available/snrna-pipeline
sudo ln -s /etc/nginx/sites-available/snrna-pipeline /etc/nginx/sites-enabled/
sudo systemctl restart nginx
```

## Configuration

### Environment Variables
```bash
# Host configuration
export SHINY_HOST=0.0.0.0
export SHINY_PORT=3838

# Memory configuration
export R_MAX_MEM_SIZE=8G

# Library path
export R_LIBS_USER=~/.local/lib/R/library
```

### Configuration Files
- `config/settings.yaml`: Pipeline configuration
- `nginx.conf`: NGINX reverse proxy configuration
- `docker-compose.yml`: Docker services configuration
- `Dockerfile`: Container build configuration

## Usage Examples

### Web Interface
1. Access the application at `http://your-server:3838`
2. Upload H5 or RDS files
3. Configure analysis parameters
4. Run the pipeline
5. View and download results

### Command Line Interface
```bash
# Basic usage with H5 file
Rscript run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name MyProject

# Advanced usage with custom parameters
Rscript run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name MyProject \
  --min_features 300 \
  --n_variable_features 3000 \
  --clustering_resolution 0.8 \
  --find_markers \
  --verbose

# Using RDS file
Rscript run_pipeline_terminal.R \
  --rds_input data/sample.rds \
  --project_name MyProject
```

### Docker Commands
```bash
# View logs
./deploy.sh logs

# Check status
./deploy.sh status

# Stop services
./deploy.sh stop

# Restart services
./deploy.sh restart

# Clean up
./deploy.sh clean
```

## Scaling Features

### 1. Server-Ready Configuration
- Uses `0.0.0.0` host instead of localhost
- Environment variable configuration
- No browser auto-launch on servers

### 2. Containerization
- Docker support for consistent deployment
- Docker Compose for multi-service deployment
- NGINX reverse proxy for production

### 3. Multiple Interfaces
- Web-based Shiny interface
- Command-line terminal interface
- Original pipeline script compatibility

### 4. Resource Management
- Configurable memory limits
- Parallel processing support
- Optimized for large datasets

### 5. Production Features
- Health checks
- Logging and monitoring
- SSL/TLS support (configured in nginx.conf)
- Rate limiting
- Gzip compression

## Troubleshooting

### Common Issues

#### Port Already in Use
```bash
# Use different port
./deploy.sh deploy -p 8080
```

#### Memory Issues
```bash
# Increase memory limit
export R_MAX_MEM_SIZE=16G
./deploy.sh deploy
```

#### Permission Errors
```bash
# Fix permissions
sudo chown -R $USER:$USER .
chmod +x deploy.sh
```

#### Docker Issues
```bash
# Clean and redeploy
./deploy.sh clean
./deploy.sh deploy
```

### Logs and Debugging
```bash
# View application logs
./deploy.sh logs

# Check service status
./deploy.sh status

# View Docker logs
docker-compose logs -f snrna-pipeline
```

## Security Considerations

### Production Deployment
1. **SSL/TLS**: Configure certificates in nginx.conf
2. **Firewall**: Open only necessary ports (80, 443, 3838)
3. **Authentication**: Implement user authentication if needed
4. **Updates**: Regularly update dependencies and system packages

### Data Security
1. **Input Validation**: All inputs are validated
2. **File Permissions**: Proper file permissions configured
3. **Temporary Files**: Cleaned up automatically
4. **Logging**: Secure logging without sensitive data

## Performance Optimization

### For Large Datasets
1. **Memory**: Increase R_MAX_MEM_SIZE
2. **CPU**: Use parallel processing
3. **Storage**: Ensure sufficient disk space
4. **Network**: Use local storage for data

### For High Traffic
1. **Load Balancing**: Use multiple instances
2. **Caching**: Implement result caching
3. **CDN**: Use content delivery network for static files
4. **Monitoring**: Monitor resource usage

## Monitoring and Maintenance

### Health Checks
```bash
# Check application health
curl http://localhost:3838

# Check Docker container health
docker ps
```

### Updates
```bash
# Update application
git pull origin main
./deploy.sh clean
./deploy.sh deploy
```

### Backup
```bash
# Backup configuration
tar -czf backup-$(date +%Y%m%d).tar.gz config/ data/ results/
```

## Support

### Documentation
- `README.md`: Main documentation
- `shiny_app/docs/`: Application documentation
- `DEPLOYMENT_GUIDE.md`: This guide

### Getting Help
1. Check the troubleshooting section
2. Review logs using `./deploy.sh logs`
3. Open an issue on GitHub
4. Contact the development team

## Conclusion

The snRNA-seq Pipeline has been successfully scaled and optimized for server deployment. Key improvements include:

1. **Fixed all Shiny app errors**
2. **Added server-ready configuration**
3. **Implemented Docker containerization**
4. **Created comprehensive terminal interface**
5. **Added production deployment features**
6. **Provided extensive documentation**

The pipeline now supports:
- Local development
- Server deployment
- Containerized deployment
- Production scaling
- Multiple user interfaces
- Comprehensive error handling

All components are tested and ready for production use.

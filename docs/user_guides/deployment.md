# Deployment Guide

This guide covers deploying the snRNA-seq pipeline in various research computing environments, from local development to production servers.

## Local Development Setup

### Prerequisites

- **R** (version 4.0 or higher)
- **RStudio** (recommended)
- **Git**
- **System dependencies** (Ubuntu/Debian):
  ```bash
  sudo apt update
  sudo apt install r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev
  ```

### Quick Setup

```bash
# Clone repository
git clone <repository-url>
cd snRNA-seq-Pipeline

# Run automated setup
./launch.sh setup

# Verify installation
Rscript scripts/run_pipeline_terminal.R --help
```

### Development Mode

```bash
# Start Shiny app for development
./launch.sh shiny

# Run terminal pipeline
./launch.sh terminal
```

## Research Computing Cluster Deployment

### SLURM Cluster Setup

Create a SLURM submission script:

```bash
#!/bin/bash
# submit_pipeline.sh

#SBATCH --job-name=snrna_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

# Load modules (adjust for your cluster)
module load R/4.2.0
module load hdf5/1.12.0

# Set environment variables
export R_MAX_MEM_SIZE=32G
export OMP_NUM_THREADS=8

# Run pipeline
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name ClusterAnalysis \
  --min_features 200 \
  --n_variable_features 2000 \
  --clustering_resolution 0.5 \
  --find_markers \
  --verbose
```

### PBS/Torque Cluster Setup

```bash
#!/bin/bash
# submit_pipeline.pbs

#PBS -N snrna_pipeline
#PBS -o logs/pipeline_$PBS_JOBID.out
#PBS -e logs/pipeline_$PBS_JOBID.err
#PBS -l walltime=24:00:00
#PBS -l mem=32gb
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR

# Load modules
module load R/4.2.0

# Run pipeline
Rscript scripts/run_pipeline_terminal.R \
  --h5_input data/sample.h5 \
  --project_name ClusterAnalysis
```

### Multi-Sample Cluster Processing

```bash
#!/bin/bash
# submit_multi_sample.sh

#SBATCH --job-name=multi_sample
#SBATCH --output=logs/multi_%j.out
#SBATCH --error=logs/multi_%j.err
#PBS -l walltime=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

# Run multi-sample analysis
Rscript scripts/run_multi_sample_pipeline.R \
  --h5_inputs data/sample1.h5 data/sample2.h5 data/sample3.h5 \
  --project_name MultiSampleCluster \
  --parallel \
  --n_cores 16 \
  --find_markers
```

## Docker Deployment

### Development Container

```bash
# Build development container
docker build -t snrna-pipeline:dev .

# Run with data volume
docker run -it --rm \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/results:/app/results \
  -p 3838:3838 \
  snrna-pipeline:dev
```

### Production Container

```bash
# Build production container
docker build -t snrna-pipeline:prod --target production .

# Run production container
docker run -d \
  --name snrna-pipeline \
  -v /path/to/data:/app/data \
  -v /path/to/results:/app/results \
  -p 3838:3838 \
  --restart unless-stopped \
  snrna-pipeline:prod
```

### Docker Compose Setup

```yaml
# docker-compose.yml
version: '3.8'

services:
  snrna-pipeline:
    build: .
    ports:
      - "3838:3838"
    volumes:
      - ./data:/app/data
      - ./results:/app/results
      - ./logs:/app/logs
    environment:
      - R_MAX_MEM_SIZE=16G
      - SHINY_HOST=0.0.0.0
      - SHINY_PORT=3838
    restart: unless-stopped

  nginx:
    image: nginx:alpine
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf
      - ./ssl:/etc/nginx/ssl
    depends_on:
      - snrna-pipeline
    restart: unless-stopped
```

## Server Deployment

### Ubuntu/Debian Server Setup

```bash
#!/bin/bash
# server_setup.sh

# Update system
sudo apt update && sudo apt upgrade -y

# Install R and dependencies
sudo apt install -y r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev

# Install Shiny Server
sudo apt install -y gdebi-core
wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb
sudo gdebi shiny-server-1.5.20.1002-amd64.deb

# Install NGINX
sudo apt install -y nginx

# Setup firewall
sudo ufw allow ssh
sudo ufw allow 'Nginx Full'
sudo ufw enable
```

### Shiny Server Configuration

```bash
# /etc/shiny-server/shiny-server.conf
server {
  listen 3838;
  
  location / {
    site_dir /srv/shiny-server/snrna-pipeline;
    log_dir /var/log/shiny-server;
    directory_index on;
  }
  
  location /snrna-pipeline {
    site_dir /srv/shiny-server/snrna-pipeline;
    log_dir /var/log/shiny-server;
  }
}
```

### NGINX Reverse Proxy

```nginx
# /etc/nginx/sites-available/snrna-pipeline
server {
    listen 80;
    server_name your-domain.com;
    
    location / {
        proxy_pass http://localhost:3838;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_cache_bypass $http_upgrade;
    }
}
```

### SSL Configuration

```bash
# Install Certbot
sudo apt install -y certbot python3-certbot-nginx

# Obtain SSL certificate
sudo certbot --nginx -d your-domain.com

# Auto-renewal
sudo crontab -e
# Add: 0 12 * * * /usr/bin/certbot renew --quiet
```

## Cloud Deployment

### AWS EC2 Deployment

```bash
#!/bin/bash
# aws_setup.sh

# Update system
sudo yum update -y

# Install R
sudo yum install -y R

# Install additional dependencies
sudo yum install -y openssl-devel libcurl-devel libxml2-devel

# Install Shiny Server
sudo yum install -y wget
wget https://download3.rstudio.org/centos7/x86_64/shiny-server-1.5.20.1002-x86_64.rpm
sudo yum install -y shiny-server-1.5.20.1002-x86_64.rpm

# Configure security groups
# Allow ports 22 (SSH), 80 (HTTP), 443 (HTTPS), 3838 (Shiny)
```

### Google Cloud Platform

```bash
#!/bin/bash
# gcp_setup.sh

# Update system
sudo apt update && sudo apt upgrade -y

# Install R
sudo apt install -y r-base r-base-dev

# Install Shiny Server
sudo apt install -y gdebi-core
wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb
sudo gdebi shiny-server-1.5.20.1002-amd64.deb

# Configure firewall
gcloud compute firewall-rules create allow-shiny \
  --allow tcp:3838 \
  --source-ranges 0.0.0.0/0 \
  --description "Allow Shiny Server"
```

### Azure VM Deployment

```bash
#!/bin/bash
# azure_setup.sh

# Update system
sudo apt update && sudo apt upgrade -y

# Install R
sudo apt install -y r-base r-base-dev

# Install Shiny Server
sudo apt install -y gdebi-core
wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb
sudo gdebi shiny-server-1.5.20.1002-amd64.deb

# Configure network security group
# Allow ports 22, 80, 443, 3838
```

## High-Performance Computing (HPC)

### Module System Configuration

Create module files for your HPC system:

```bash
# /usr/local/modulefiles/R/4.2.0
#%Module1.0
##
## R 4.2.0 modulefile
##

proc ModulesHelp { } {
    puts stderr "This module sets up the environment for R 4.2.0"
}

module-whatis "Sets up the environment for R 4.2.0"

prepend-path    PATH            /usr/local/R-4.2.0/bin
prepend-path    LD_LIBRARY_PATH /usr/local/R-4.2.0/lib64/R/lib
prepend-path    MANPATH         /usr/local/R-4.2.0/share/man
```

### Resource Management

```bash
#!/bin/bash
# resource_monitor.sh

# Monitor resource usage during pipeline execution
while true; do
    echo "$(date): Memory usage: $(free -h | grep Mem | awk '{print $3"/"$2}')"
    echo "$(date): CPU usage: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}')"
    echo "$(date): Disk usage: $(df -h / | tail -1 | awk '{print $5}')"
    sleep 60
done
```

## Monitoring and Logging

### Log Configuration

```r
# logging_config.R
library(logger)

# Configure logging
log_threshold(INFO)
log_appender(appender_file("logs/pipeline.log"))

# Log pipeline events
log_info("Pipeline started")
log_info("Processing sample: %s", sample_name)
log_info("Pipeline completed successfully")
```

### System Monitoring

```bash
#!/bin/bash
# monitor_system.sh

# Monitor system resources
while true; do
    echo "=== $(date) ===" >> logs/system_monitor.log
    echo "Memory: $(free -h)" >> logs/system_monitor.log
    echo "CPU: $(top -bn1 | grep "Cpu(s)")" >> logs/system_monitor.log
    echo "Disk: $(df -h)" >> logs/system_monitor.log
    echo "Processes: $(ps aux | grep R | wc -l)" >> logs/system_monitor.log
    echo "" >> logs/system_monitor.log
    sleep 300
done
```

## Backup and Recovery

### Data Backup Strategy

```bash
#!/bin/bash
# backup_data.sh

# Backup configuration and results
tar -czf backups/config_$(date +%Y%m%d).tar.gz config/
tar -czf backups/results_$(date +%Y%m%d).tar.gz results/

# Sync to remote storage
rsync -av backups/ user@backup-server:/backups/snrna-pipeline/

# Clean old backups (keep 30 days)
find backups/ -name "*.tar.gz" -mtime +30 -delete
```

### Recovery Procedures

```bash
#!/bin/bash
# recovery.sh

# Restore from backup
tar -xzf backups/config_$(date +%Y%m%d).tar.gz
tar -xzf backups/results_$(date +%Y%m%d).tar.gz

# Restart services
sudo systemctl restart shiny-server
sudo systemctl restart nginx
```

## Security Considerations

### User Access Control

```bash
# Create dedicated user
sudo useradd -m -s /bin/bash snrna
sudo usermod -aG sudo snrna

# Set up SSH keys
sudo mkdir -p /home/snrna/.ssh
sudo cp ~/.ssh/authorized_keys /home/snrna/.ssh/
sudo chown -R snrna:snrna /home/snrna/.ssh
```

### Firewall Configuration

```bash
# Configure UFW firewall
sudo ufw default deny incoming
sudo ufw default allow outgoing
sudo ufw allow ssh
sudo ufw allow 'Nginx Full'
sudo ufw allow 3838
sudo ufw enable
```

### SSL/TLS Configuration

```nginx
# Enhanced SSL configuration
ssl_protocols TLSv1.2 TLSv1.3;
ssl_ciphers ECDHE-RSA-AES256-GCM-SHA512:DHE-RSA-AES256-GCM-SHA512:ECDHE-RSA-AES256-GCM-SHA384:DHE-RSA-AES256-GCM-SHA384;
ssl_prefer_server_ciphers off;
ssl_session_cache shared:SSL:10m;
ssl_session_timeout 10m;
```

## Performance Tuning

### R Performance Optimization

```r
# R performance settings
options(mc.cores = parallel::detectCores())
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# Memory management
gc()
memory.limit(size = 32000)  # 32GB
```

### System Performance Tuning

```bash
# Kernel parameters for high memory usage
echo 'vm.max_map_count=262144' | sudo tee -a /etc/sysctl.conf
echo 'vm.swappiness=10' | sudo tee -a /etc/sysctl.conf
sudo sysctl -p
```

This deployment guide provides comprehensive instructions for various research computing environments. Choose the appropriate deployment method based on your infrastructure and requirements.

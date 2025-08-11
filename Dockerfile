# Use R base image with Ubuntu
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV R_VERSION=4.3.2
ENV SHINY_HOST=0.0.0.0
ENV SHINY_PORT=3838

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    gnupg2 \
    software-properties-common \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libgsl-dev \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libudunits2-dev \
    libgmp3-dev \
    libmpfr-dev \
    libfftw3-dev \
    libgsl-dev \
    libhdf5-dev \
    libnetcdf-dev \
    libopenblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    gfortran \
    git \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Add CRAN repository
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | apt-key add - \
    && echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" >> /etc/apt/sources.list

# Install R
RUN apt-get update && apt-get install -y \
    r-base=${R_VERSION}-1ubuntu2 \
    r-base-dev=${R_VERSION}-1ubuntu2 \
    && rm -rf /var/lib/apt/lists/*

# Install Shiny Server
RUN wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb \
    && apt-get update \
    && apt-get install -y gdebi-core \
    && gdebi --non-interactive shiny-server-1.5.20.1002-amd64.deb \
    && rm shiny-server-1.5.20.1002-amd64.deb \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /app

# Copy R scripts and configuration
COPY . .

# Install R packages
RUN Rscript setup/install_dependencies.R

# Install additional packages for multi-sample processing and Shiny
RUN R -e "install.packages(c('future', 'future.apply', 'parallel', 'purrr', 'shinyFiles', 'shinyWidgets', 'shinyjs', 'DT', 'fs', 'promises', 'Matrix', 'hdf5r', 'stringr'), repos='https://cran.rstudio.com/')"

# Make scripts executable
RUN chmod +x run_shiny_app.R \
    && chmod +x run_pipeline_terminal.R \
    && chmod +x seurat_pipeline.R

# Create necessary directories
RUN mkdir -p /var/log/shiny-server \
    && mkdir -p /var/lib/shiny-server \
    && mkdir -p /srv/shiny-server \
    && mkdir -p /app/logs \
    && mkdir -p /app/results \
    && mkdir -p /app/temp

# Copy app to Shiny Server directory
RUN cp -r shiny_app /srv/shiny-server/snrna-pipeline

# Configure Shiny Server
RUN echo 'server {\n\
    listen 3838;\n\
    location / {\n\
        site_dir /srv/shiny-server;\n\
        log_dir /var/log/shiny-server;\n\
        directory_index on;\n\
    }\n\
}' > /etc/shiny-server/shiny-server.conf

# Expose port
EXPOSE 3838

# Create startup script
RUN echo '#!/bin/bash\n\
echo "Starting snRNA-seq Pipeline..."\n\
echo "Available options:"\n\
echo "1. Run Shiny app: Rscript run_shiny_app.R"\n\
echo "2. Run terminal pipeline: Rscript run_pipeline_terminal.R --help"\n\
echo "3. Run original pipeline: Rscript seurat_pipeline.R --help"\n\
echo ""\n\
echo "Starting Shiny Server..."\n\
exec shiny-server' > /app/start.sh \
    && chmod +x /app/start.sh

# Set working directory
WORKDIR /app

# Default command
CMD ["/app/start.sh"]

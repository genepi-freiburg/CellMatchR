FROM rocker/shiny:4.4

# System dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN R -e "install.packages(c( \
    'shinydashboard', 'DT', 'ggplot2', 'ggfortify', 'ggrepel', \
    'ggpubr', 'ggdendro', 'gridExtra', 'readxl', 'shinyWidgets', \
    'shinyjs', 'circlize', 'viridis', 'dplyr', 'openxlsx', 'writexl' \
), repos='https://cloud.r-project.org')"

# Install Bioconductor package
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    BiocManager::install('ComplexHeatmap', ask=FALSE, update=FALSE)"

# Copy app
RUN rm -rf /srv/shiny-server/*
COPY . /srv/shiny-server/

# Fix permissions
RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 3841

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=3841)"]

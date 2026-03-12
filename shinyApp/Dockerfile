FROM rocker/shiny:latest

ENV MAKEFLAGS="-j4"

# System dependencies for R packages (one per line for better caching)
RUN apt-get update && apt-get install -y cmake && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libcurl4-openssl-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libssl-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libxml2-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libfontconfig1-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libfreetype6-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libharfbuzz-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libfribidi-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libpng-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libtiff5-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y libjpeg-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y zlib1g-dev && rm -rf /var/lib/apt/lists/*

# Install CRAN packages (one per line for better caching)
RUN R -e "install.packages('shinydashboard', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('DT', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('ggplot2', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('ggfortify', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('ggrepel', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('ggpubr', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('ggdendro', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('gridExtra', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('readxl', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('shinyWidgets', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('shinyjs', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('circlize', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('viridis', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('dplyr', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('openxlsx', repos='https://cloud.r-project.org', Ncpus=4)"
RUN R -e "install.packages('writexl', repos='https://cloud.r-project.org', Ncpus=4)"

# Install Bioconductor package
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    BiocManager::install('ComplexHeatmap', ask=FALSE, update=FALSE, Ncpus=4)"

# Copy app
RUN rm -rf /srv/shiny-server/*
COPY . /srv/shiny-server/

# Fix permissions
RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 3841

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=3841)"]

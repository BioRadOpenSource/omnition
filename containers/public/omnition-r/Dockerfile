# syntax=docker/dockerfile:1.2

# Copyright 2020 Bio-Rad Laboratories, Inc.

# Set up image version arg and cont dir arg
ARG IMAGE_VERSION
ARG CONT_DIR
ARG STAGE=""

# Starting with image already containing mamba
FROM bioraddbg/omnition-base${STAGE}:${IMAGE_VERSION}

# Building image using bash shell
SHELL ["/bin/bash", "-c"]

# Name of mamba environment being created (must match env.yaml filename)
ARG MAMBA_ENV="omnition-r"

# Building the mamba environment
USER root
COPY $CONT_DIR/$MAMBA_ENV.yaml /opt/biorad/env/
RUN micromamba create -f /opt/biorad/env/$MAMBA_ENV.yaml \
    && micromamba clean -afy

# Removing default activation commands and setting correct ones
RUN rm /root/.bashrc
RUN echo "source /etc/container.bashrc" >> /etc/bash.bashrc \
    && echo "set +u" > /etc/container.bashrc \
    && echo 'eval "$(micromamba shell hook --shell=bash)"' >> /etc/container.bashrc \
    && echo "micromamba activate $MAMBA_ENV" >> /etc/container.bashrc

# Activating environment when using non-login, non-interactive shell
ENV BASH_ENV /etc/container.bashrc
ENV ENV /etc/container.bashrc

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

# ============================================================================
# FIX FOR APPTAINER: Redirect R temp/cache directories to /tmp
# This prevents issues when home directory is not mounted in Apptainer
# Fixes segfaults and permission errors in R scripts (e.g., rnaCombineReadCounts.R)
# ============================================================================

# Set R environment variables to use /tmp for cache and temporary files
ENV R_USER_CACHE_DIR=/tmp/Rcache
ENV R_USER_DATA_DIR=/tmp/Rdata
ENV TMPDIR=/tmp
ENV TEMP=/tmp
ENV TMP=/tmp

# Create the cache directories and make them writable
RUN mkdir -p ${R_USER_CACHE_DIR} ${R_USER_DATA_DIR} && chmod -R 775 ${R_USER_CACHE_DIR} ${R_USER_DATA_DIR}

# Set R library path to include a writable location in /tmp as fallback
# NOTE: Don't set R_LIBS_USER to /tmp during build - packages would be lost!
# This will be set at runtime for Apptainer compatibility
RUN mkdir -p /tmp/Rlib && chmod -R 775 /tmp/Rlib

# Create R profile that sets these options at R startup
# Find R installation and create Rprofile.site in the correct location
RUN R_HOME=$(R RHOME) && \
    mkdir -p ${R_HOME}/etc && \
    echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' > ${R_HOME}/etc/Rprofile.site && \
    echo 'Sys.setenv(TMPDIR="/tmp")' >> ${R_HOME}/etc/Rprofile.site && \
    echo 'Sys.setenv(R_USER_CACHE_DIR="/tmp/Rcache")' >> ${R_HOME}/etc/Rprofile.site && \
    echo 'Sys.setenv(R_USER_DATA_DIR="/tmp/Rdata")' >> ${R_HOME}/etc/Rprofile.site

# ============================================================================

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/

# Install GitHub packages within the activated conda environment
# Explicitly activate the environment and install to the conda R library
RUN eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate omnition-r && \
    R -e ".libPaths(); devtools::install_github('immunogenomics/presto@31dc97fed5e2e7fc323ae4af62f72181cc51d9a3', upgrade='never', repos=NULL, lib=.libPaths()[1])" && \
    R -e ".libPaths(); devtools::install_github('satijalab/seurat-wrappers@d28512f804d5fe05e6d68900ca9221020d52cf1d', upgrade='never', repos=NULL, lib=.libPaths()[1])" && \
    R -e "cat('Library paths:\n'); print(.libPaths()); cat('\nVerifying installations:\n'); library(presto); library(SeuratWrappers); cat('SUCCESS: All packages loaded\n')"

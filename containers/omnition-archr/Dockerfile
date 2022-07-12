# syntax=docker/dockerfile:1.2

# Copyright 2020 Bio-Rad Laboratories, Inc.

# Starting with image already containing conda
FROM continuumio/miniconda3:4.7.12

# Building image using bash shell
SHELL ["/bin/bash", "-c"]

# Name of conda environment being created (must match env.yaml filename)
ARG CONDA_ENV="omnition-archr"

# Installing image dependencies
RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y procps=2:3.3.15-2 && \
    apt-get clean -y

# Building the conda environment
COPY $CONDA_ENV.yaml /opt/biorad/env/
RUN conda env create -f /opt/biorad/env/$CONDA_ENV.yaml && \
    conda clean -afy

# Removing default acitvation commands and setting correct ones
RUN rm root/.bashrc
RUN echo "source /etc/container.bashrc" >> /etc/bash.bashrc && \
    echo "set +u" > /etc/container.bashrc && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/container.bashrc && \
    echo "conda activate $CONDA_ENV" >> /etc/container.bashrc

# Adding in an empty folder to the BSgenome package which is required but not always installed depending on the source
RUN mkdir -p /opt/conda/envs/${CONDA_ENV}/lib/R/library/BSgenome/pkgtemplates/BSgenome_datapkg/inst/extdata

# Activating environment when using non-login, non-interactive shell
ENV BASH_ENV /etc/container.bashrc
ENV ENV /etc/container.bashrc

# Install Archr prereqs
RUN apt-get install -y libcairo2-dev && \
    Rscript -e 'BiocManager::install("TFBSTools")' && \
    Rscript -e 'devtools::install_github("immunogenomics/presto", ref = "e7c9c6b708584aa23e590c5722bb5d35442f183d", upgrade = "never" , repos = BiocManager::repositories())' && \
    apt-get install -y libxt-dev &&\
    apt-get clean -y

# Install ArchR
RUN Rscript -e 'devtools::install_github("GreenleafLab/ArchR", ref = "48dccf361bab6d9191a4d9b61e1317acfa71167c", upgrade = "never", repos = BiocManager::repositories())'
# RUN Rscript -e 'ArchR::installExtraPackages()'

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/


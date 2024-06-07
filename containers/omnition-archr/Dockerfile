# syntax=docker/dockerfile:1.2

# Copyright 2020 Bio-Rad Laboratories, Inc.

# Starting with image already containing mamba
FROM mambaorg/micromamba:1.0-bullseye-slim

# Get container directory as an ARG
ARG CONT_DIR

# Building image using bash shell
SHELL ["/bin/bash", "-c"]

# Name of mamba environment being created (must match env.yaml filename)
ARG MAMBA_ENV="omnition-archr"

# Installing image dependencies
USER root
RUN apt-get --allow-releaseinfo-change update \
    && apt-get install -y procps=2:3.3.17-5 gawk=1:5.1.0-1 \
    && apt-get clean -y

# Building the mamba environment
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

# Adding empty folder required for the BSgenome package
RUN mkdir -p /opt/conda/envs/${MAMBA_ENV}/lib/R/library/BSgenome/pkgtemplates/BSgenome_datapkg/inst/extdata

# Install ArchR
RUN Rscript -e 'devtools::install_github("GreenleafLab/ArchR", ref = "v1.0.2", upgrade = "never", repos = BiocManager::repositories())'
RUN apt-get update \
    && apt-get install -y libcairo2-dev \
    && Rscript -e 'install.packages("Cairo", repos = "http://cran.us.r-project.org")' \
    && Rscript -e 'ArchR::installExtraPackages()'

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/

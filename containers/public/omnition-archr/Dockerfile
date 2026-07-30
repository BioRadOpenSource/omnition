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
ARG MAMBA_ENV="omnition-archr"

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

# Adding empty folder required for the BSgenome package
RUN mkdir -p /opt/conda/envs/${MAMBA_ENV}/lib/R/library/BSgenome/pkgtemplates/BSgenome_datapkg/inst/extdata

# Install TFMPvalue from GitHub (no longer available in Bioconductor)
# Note: Modifying R version requirement for compatibility with R 4.2.2
RUN eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate ${MAMBA_ENV} && \
    cd /tmp && \
    git clone https://github.com/ge11232002/TFMPvalue.git && \
    cd TFMPvalue && \
    git checkout ddc1173 && \
    sed -i 's/R (>= 4.5.0)/R (>= 4.2.0)/g' DESCRIPTION && \
    Rscript -e 'devtools::install(".", upgrade = "never")' && \
    cd / && rm -rf /tmp/TFMPvalue

# Install ArchR (requires conda environment activation)
# Note: Installing GlobalOptions v0.1.2 for compatibility with R 4.2.2
RUN eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate ${MAMBA_ENV} && \
    Rscript -e 'devtools::install_version("GlobalOptions", version = "0.1.2", repos = "https://cloud.r-project.org", upgrade = "never")' && \
    Rscript -e 'devtools::install_github("GreenleafLab/ArchR", ref = "v1.0.2", upgrade = "never", repos = BiocManager::repositories())'

# Install Cairo and ArchR extra packages (requires conda environment activation)
RUN apt-get update && \
    apt-get install -y libcairo2-dev && \
    eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate ${MAMBA_ENV} && \
    Rscript -e 'install.packages("Cairo", repos = "http://cran.us.r-project.org")' && \
    Rscript -e 'library(ArchR); ArchR::installExtraPackages()'

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/

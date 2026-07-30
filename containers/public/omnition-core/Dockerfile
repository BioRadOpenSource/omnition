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
ARG MAMBA_ENV="omnition-core"

# Building the mamba environment
# Note: This build must occur from the repository root
USER root
COPY $CONT_DIR/$MAMBA_ENV.yaml /opt/biorad/env/
RUN micromamba create -f /opt/biorad/env/$MAMBA_ENV.yaml \
    && micromamba clean -afy

# Install bragg
WORKDIR /
RUN --mount=type=secret,id=GITHUB_TOKEN git clone --branch v0.4.2 \
    https://$(cat /run/secrets/GITHUB_TOKEN):x-oauth-basic@github.com/digitalbiology/bragg.git \
    && cd ./bragg \
    && eval "$(micromamba shell hook --shell=bash)" \
    && micromamba activate $MAMBA_ENV \
    && pip install -e ./

# Removing default activation commands and setting correct ones
RUN rm /root/.bashrc
RUN echo "source /etc/container.bashrc" >> /etc/bash.bashrc \
    && echo "set +u" > /etc/container.bashrc \
    && echo 'eval "$(micromamba shell hook --shell=bash)"' >> /etc/container.bashrc \
    && echo "micromamba activate $MAMBA_ENV" >> /etc/container.bashrc

# Activating environment when using non-login, non-interactive shell
ENV BASH_ENV /etc/container.bashrc
ENV ENV /etc/container.bashrc

# Copy assets to container
# Note: This build must occur from the repository root
COPY $CONT_DIR/assets /opt/biorad/assets/
RUN ls --recursive /opt/biorad/

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/

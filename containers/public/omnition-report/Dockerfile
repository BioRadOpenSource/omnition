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
ARG MAMBA_ENV="omnition-report"

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

# ============================================================================
# FIX FOR APPTAINER: Redirect cache directories to /tmp
# This prevents npm and matplotlib errors when home directory is not mounted
# Fixes EROFS errors and matplotlib config warnings
# ============================================================================

# Set npm to use /tmp for cache instead of ~/.npm
ENV NPM_CONFIG_CACHE=/tmp/npm_cache
ENV NPM_CONFIG_PREFIX=/tmp/npm_global
ENV NPM_CONFIG_TMP=/tmp/npm_tmp

# Set matplotlib to use /tmp for config instead of ~/.config/matplotlib
ENV MPLCONFIGDIR=/tmp/matplotlib

# General temp directory settings
ENV TMPDIR=/tmp
ENV TEMP=/tmp
ENV TMP=/tmp

# Set XDG directories to /tmp to catch other cache/config attempts
ENV XDG_CACHE_HOME=/tmp/cache
ENV XDG_CONFIG_HOME=/tmp/config
ENV XDG_DATA_HOME=/tmp/data

# Create cache directories with proper permissions
RUN mkdir -p ${NPM_CONFIG_CACHE} ${NPM_CONFIG_PREFIX} ${NPM_CONFIG_TMP} ${MPLCONFIGDIR} && \
    chmod -R 775 ${NPM_CONFIG_CACHE} ${NPM_CONFIG_PREFIX} ${NPM_CONFIG_TMP} ${MPLCONFIGDIR}

RUN mkdir -p ${XDG_CACHE_HOME} ${XDG_CONFIG_HOME} ${XDG_DATA_HOME} && \
    chmod -R 775 ${XDG_CACHE_HOME} ${XDG_CONFIG_HOME} ${XDG_DATA_HOME}

# ============================================================================

# Install report application
ENV APP_HOME /opt/biorad/app/
WORKDIR $APP_HOME

# Install node prereqs
# Ref: https://deb.nodesource.com/setup_16.x
RUN curl -sL https://deb.nodesource.com/setup_16.x | bash -

# Install application
RUN --mount=type=secret,id=GITHUB_TOKEN micromamba activate $MAMBA_ENV \
    && git clone --branch v1.1.0-beta.6 \
    https://$(cat /run/secrets/GITHUB_TOKEN):x-oauth-basic@github.com/digitalbiology/single-cell-report.git ./ \
    && for FILE in $(cat .dockerignore); do rm -rf $FILE; done \
    && rm -rf .dockerignore .git* \
    && npm i --cache ${NPM_CONFIG_CACHE}

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/

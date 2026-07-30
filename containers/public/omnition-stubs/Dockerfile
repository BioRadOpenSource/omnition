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

# Note: This build must occur from the repository root
# Copying stub_functions.sh into the image
USER root
COPY $CONT_DIR/assets /opt/biorad/assets/

# Removing default activation commands and setting correct ones
RUN rm /root/.bashrc
RUN echo "source /etc/container.bashrc" >> /etc/bash.bashrc \
    && echo "set +u" > /etc/container.bashrc \
    && echo "source /opt/biorad/assets/stub_functions.sh" >> /etc/container.bashrc

# Activating environment when using non-login, non-interactive shell
ENV BASH_ENV /etc/container.bashrc
ENV ENV /etc/container.bashrc

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/

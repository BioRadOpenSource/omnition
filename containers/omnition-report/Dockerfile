# syntax=docker/dockerfile:1.2

# Copyright 2020 Bio-Rad Laboratories, Inc.

# Starting with image already containing mamba
FROM mambaorg/micromamba:1.0-bullseye-slim

# Building image using bash shell
SHELL ["/bin/bash", "-c"]

# Name of mamba environment being created (must match env.yaml filename)
ARG MAMBA_ENV="omnition-report"

# Installing image dependencies
USER root
RUN apt-get --allow-releaseinfo-change update \
    && apt-get install -y procps=2:3.3.17-5 gawk=1:5.1.0-1 \
    && apt-get clean -y

# Building the mamba environment
COPY $MAMBA_ENV.yaml /opt/biorad/env/
RUN micromamba create -f /opt/biorad/env/$MAMBA_ENV.yaml \
    && micromamba clean -afy

# Removing default acitvation commands and setting correct ones
RUN rm /root/.bashrc
RUN echo "source /etc/container.bashrc" >> /etc/bash.bashrc \
    && echo "set +u" > /etc/container.bashrc \
    && echo 'eval "$(micromamba shell hook --shell=bash)"' >> /etc/container.bashrc \
    && echo "micromamba activate $MAMBA_ENV" >> /etc/container.bashrc

# Activating environment when using non-login, non-interactive shell
ENV BASH_ENV /etc/container.bashrc
ENV ENV /etc/container.bashrc

# Install report application
ENV APP_HOME /opt/biorad/app/
WORKDIR $APP_HOME 

# Install node prereqs
# Ref: https://deb.nodesource.com/setup_12.x
RUN curl -sL https://deb.nodesource.com/setup_12.x | bash -

# Install application
RUN --mount=type=secret,id=GITHUB_TOKEN micromamba activate $MAMBA_ENV \
    && git clone --branch v1.0.3-alpha.1 \
    https://$(cat /run/secrets/GITHUB_TOKEN):x-oauth-basic@github.com/BioRadOpenSource/cuttlesoft-single-cell-report.git ./ \
    && for FILE in $(cat .dockerignore); do rm -rf $FILE; done \
    && rm -rf .dockerignore .git* \
    && npm i

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/

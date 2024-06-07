# syntax=docker/dockerfile:1.2

# Copyright 2020 Bio-Rad Laboratories, Inc.

FROM rust:1.69 as dead-build

# Build DEAD dependencies to cache them
WORKDIR /
RUN USER=root cargo new --bin debarcoder

# Install DEAD
WORKDIR /dead
RUN --mount=type=secret,id=GITHUB_TOKEN git clone --branch v0.2.22 \
    https://$(cat /run/secrets/GITHUB_TOKEN):x-oauth-basic@github.com/digitalbiology/debarcoder.git \
    && cd ./debarcoder \
    && cargo build --release \
    && rm src/*.rs  

# Copyright 2020 Bio-Rad Laboratories, Inc.
FROM rust:1.70 as omnition-build

# Build DEAD dependencies to cache them
RUN apt-get update \
    && apt-get install -y gcc-multilib \
    && apt-get install -y g++ cmake libsasl2-dev build-essential \
    && apt-get install -y linux-libc-dev zlib1g zlib1g-dev

WORKDIR /
RUN --mount=type=secret,id=GITHUB_TOKEN git clone --branch develop \
    https://$(cat /run/secrets/GITHUB_TOKEN):x-oauth-basic@github.com/digitalbiology/omnition-rs.git \
    && cd ./omnition-rs \
    && cargo build --release \
    && rm src/*.rs  

# Starting with image already containing mamba
FROM mambaorg/micromamba:1.0-bullseye-slim

# Get container directory as an ARG
ARG CONT_DIR

# Building image using bash shell
SHELL ["/bin/bash", "-c"]

# Name of mamba environment being created (must match env.yaml filename)
ARG MAMBA_ENV="omnition-dbg"

# Installing image dependencies
USER root
RUN apt-get --allow-releaseinfo-change update \
    && apt-get install -y procps=2:3.3.17-5 gawk=1:5.1.0-1 \
    && apt-get install -y openssh-client \
    && apt-get clean -y

# Building the mamba environment
COPY $CONT_DIR/$MAMBA_ENV.yaml /opt/biorad/env/
RUN micromamba create -f /opt/biorad/env/$MAMBA_ENV.yaml \
    && micromamba clean -afy

RUN rm /root/.bashrc
RUN echo "source /etc/container.bashrc" >> /etc/bash.bashrc \
    && echo "set +u" > /etc/container.bashrc \
    && echo 'eval "$(micromamba shell hook --shell=bash)"' >> /etc/container.bashrc \
    && echo "micromamba activate $MAMBA_ENV" >> /etc/container.bashrc

# Activating environment when using non-login, non-interactive shell
ENV BASH_ENV /etc/container.bashrc
ENV ENV /etc/container.bashrc

# Creating bin directory structure
RUN mkdir -p /opt/biorad/bin \
    && mkdir -p /opt/biorad/assets/dead \
    && mkdir -p /usr/local/src/dead/bin \
    && mkdir -p /usr/local/src/dead/debarcoder/configs \
    && mkdir -p /usr/local/src/omnition-rs/bin

# Copying binaries from previous stages of the build
COPY --from=dead-build /dead/debarcoder/target/release/dead /usr/local/src/dead/bin/
COPY --from=dead-build /dead/debarcoder/configs/* /usr/local/src/dead/debarcoder/configs/
COPY --from=omnition-build /omnition-rs/target/release/omnition /usr/local/src/omnition-rs/bin/

# Adding symbolic links to Bio-Rad bin
RUN ln -s /usr/local/src/dead/bin/dead /opt/biorad/bin/ \
    && ln -s /usr/local/src/omnition-rs/bin/omnition /opt/biorad/bin/ \
    && ln -s /usr/local/src/dead/debarcoder/configs/* /opt/biorad/assets/dead/

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

# Setting default working directory
# WARNING: This is not preserved by Singularity
WORKDIR /opt/biorad/

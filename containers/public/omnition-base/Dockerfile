# syntax=docker/dockerfile:1.2

# Copyright 2020 Bio-Rad Laboratories, Inc.

# Starting with image already containing mamba
FROM mambaorg/micromamba:bookworm-slim

# Building image using bash shell
SHELL ["/bin/bash", "-c"]

# Installing image dependencies
USER root
RUN apt-get --allow-releaseinfo-change update \
    && apt-get install -y procps=2:4.0.2-3 gawk=1:5.2.1-2 git=1:2.39.5-0+deb12u2 \
    && apt-get clean -y

# Adding Seqera Platform bin to path
ENV PATH /usr/local/bin/omnition:/opt/biorad/nextflow-bin/omnition:$PATH

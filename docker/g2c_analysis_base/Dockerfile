# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This Dockerfile builds the base of the image used for all analysis after initial
# variant calling (e.g., starting with quality control/filtering and proceeding
# through final formal analyses)

# Base image: conda-forge's miniconda build with mamba pre-installed (for building python environment prior to conda-pack)
FROM condaforge/mambaforge:24.9.2-0
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Install git & curl
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends git curl unzip

# Clone pancan_germline_wgs repo and check out specified hash
# Note that this is only used for pulling the version-tracked conda environment .yaml
# It will be re-cloned when building the child image, g2c_analysis
ARG G2C_TAG=main
RUN git clone \
    https://github.com/collins-genomics/pancan_germline_wgs.git \
    opt/pancan_germline_wgs && \
    cd opt/pancan_germline_wgs && \
    git pull && \
    git checkout $G2C_TAG && \
    cd -

# Create default python environment (environment must be present in build context)
ENV CONDA_SUBDIR=linux-64
RUN conda config --env --set subdir linux-64 && \
    mamba update --all && \
    mamba env create --file=/opt/pancan_germline_wgs/refs/config/environment.g2c_analysis.yml

# Manually install packages (only necessary for local debugging/environment rework)
# For whatever reason, it seems like zip and less are easier to install via conda than mamba
# conda install -y conda-forge::zip conda-forge::less
# mamba install -y bioconda::bcftools bioconda::pysam bioconda::tabix bioconda::cromshell \
#                  conda-forge::biopython conda-forge::conda-pack bioconda::pybedtools \
#                  bioconda::pybigwig
# mamba install -y conda-forge::pandas conda-forge::numpy conda-forge::scipy pip==25.0 \
#                  conda-forge::scikit-learn conda-forge::cython
# mamba install -y conda-forge::r-base conda-forge::r-desctools r::r-eql conda-forge::r-argparse \
#                 conda-forge::r-bedr conda-forge::r-beeswarm conda-forge::r-caret \
#                 conda-forge::r-readxl conda-forge::r-survival conda-forge::r-vioplot \
#                 conda-forge::r-viridis conda-forge::r-optparse
# mamba install -y conda-forge::google-cloud-sdk conda-forge::r-glmnet

# Use conda-pack to create a standalone enviornment
RUN mamba install conda-pack && \
    conda-pack -n g2c -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar -xf /tmp/env.tar && \
    rm /tmp/env.tar && \
    /venv/bin/conda-unpack

# Launch bash as entrypoint
CMD ["/bin/bash"]

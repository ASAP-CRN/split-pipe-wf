
ARG IMAGE_NAME
ARG IMAGE_TAG
ARG MINICONDA_VERSION

# Unzip the split-pipe source code and copy it into the image
FROM ubuntu:jammy AS tools

ARG SPLITPIPE_VERSION

RUN apt-get -yq update \
    && apt-get -yq install \
        unzip \
        wget \
    && rm -rf /var/lib/apt/lists/*

COPY split-pipe-v${SPLITPIPE_VERSION}.zip /tmp/split-pipe.zip
RUN unzip /tmp/split-pipe.zip -d /tmp/split-pipe-extracted \
    && mv /tmp/split-pipe-extracted/*/ /opt/split-pipe \
    && rm -rf /tmp/split-pipe-extracted

# Installing Miniconda for split-pipe conda environment
FROM anaconda/miniconda:${MINICONDA_VERSION} AS base 

LABEL MAINTAINER="Alexander Morin <alexander.morin@dnastack.com>"

ARG PYTHON_VERSION
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true

RUN apt-get -yq update \
    && apt-get -yq install gcc \
    && rm -rf /var/lib/apt/lists/*

RUN conda create -n split-pipe python=${PYTHON_VERSION} -y && \
    conda config --set default_activation_env split-pipe

# Copy and install split-pipe
COPY --from=tools /opt/split-pipe /opt/split-pipe

RUN conda run -n split-pipe bash /opt/split-pipe/install_dependencies.sh -i -y \
    && conda run -n split-pipe pip install --no-cache-dir /opt/split-pipe/ \
    && ln -s /opt/split-pipe/split-pipe /usr/local/bin/split-pipe
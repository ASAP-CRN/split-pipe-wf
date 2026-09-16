ARG MINICONDA_VERSION

# Unzip the split-pipe source code and copy it into the image
FROM ubuntu:jammy AS tools

ARG SPLIT_PIPE_VERSION

RUN apt-get -yq update \
    && apt-get -yq install unzip \
    && rm -rf /var/lib/apt/lists/*

COPY split-pipe-v${SPLIT_PIPE_VERSION}.zip /tmp/split-pipe.zip
RUN unzip /tmp/split-pipe.zip -d /tmp/split-pipe-extracted \
    && mv /tmp/split-pipe-extracted/*/ /opt/split-pipe \
    && rm -rf /tmp/split-pipe-extracted

# Installing Miniconda for split-pipe conda installation
FROM anaconda/miniconda:${MINICONDA_VERSION} AS base

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

LABEL maintainer="bioinformatics@dnastack.com"

ARG PYTHON_VERSION
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true

RUN apt-get -yq update \
    && apt-get -yq install \
        sudo \
        gcc \
        wget \
    && rm -rf /var/lib/apt/lists/*

RUN conda create -n split-pipe python=${PYTHON_VERSION} -y

# Copy and install split-pipe
COPY --from=tools /opt/split-pipe /opt/split-pipe

RUN conda run -n split-pipe bash /opt/split-pipe/install_dependencies.sh -i -y \
    && conda run -n split-pipe pip install --no-cache-dir /opt/split-pipe/ \
    && conda install -n split-pipe -y bioconda::fastqc \
    && ln -s /opt/miniconda3/envs/split-pipe/bin/split-pipe /usr/local/bin/split-pipe

# Gcloud for uploading output
ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"
RUN wget "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" \
    && tar -zxvf "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" --directory /opt \
    && rm "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz"

ENV PATH "${PATH}:/opt/google-cloud-sdk/bin"
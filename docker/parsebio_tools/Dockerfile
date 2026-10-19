ARG SC_TOOLS_DIGEST
FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/sc_tools@${SC_TOOLS_DIGEST}

ARG IMAGE_NAME
ENV IMAGE_NAME="${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG="${IMAGE_TAG}"
ARG SC_TOOLS_VERSION
ARG SC_TOOLS_DIGEST

LABEL maintainer="bioinformatics@dnastack.com"
LABEL base.sc_tools.version="${SC_TOOLS_VERSION}" \
      base.sc_tools.digest="${SC_TOOLS_DIGEST}"

# Parse Bio postprocess scripts; sc_tools' own scripts (merge_and_plot_qc, filter, ...) stay on PATH
COPY scripts /opt/scripts/parsebio
RUN chmod +x /opt/scripts/parsebio/*.py
ENV PATH="/opt/scripts/parsebio:${PATH}"

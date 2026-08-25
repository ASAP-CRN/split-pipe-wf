#!/bin/bash
set -euo pipefail

script_dir=$(dirname "$(readlink -f "$0")")

# sets split-pipe version
. "${script_dir}/build.env"

bucket_path="gs://dnastack-split-pipe/split-pipe-v${SPLITPIPE_VERSION}.zip"
dest="${script_dir}/split-pipe-v${SPLITPIPE_VERSION}.zip"

echo "Fetching split-pipe from ${bucket_path}"
gcloud storage cp "${bucket_path}" "${dest}"
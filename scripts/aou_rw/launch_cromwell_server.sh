#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Cromwell server launch script for AoU RW v2.0 (Verily Pre)

# Note that this must be run in a separate terminal

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://rw-migration-aou-rw-84a0039b

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Set terminal timezone to Eastern US time
export TZDIR=$(python -c "import tzdata, os; print(os.path.join(os.path.dirname(tzdata.__file__), 'zoneinfo'))")
export TZ=America/New_York
date +"%Y-%m-%d %H:%M:%S %Z %z"

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Reassign $WORKSPACE_BUCKET to match AoU v1-migrated URIs
case "$WN" in
  "w1")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-84a0039b"
    ;;
  "w2")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-78e2871d"
    ;;
  "w3")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-3c78b3b7"
    ;;
  "w4")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-efb2fd38"
    ;;
  "w5")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-484d2a66"
    ;;
  "dev")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-e34d8d8a"
    ;;
  *)
    echo "UNKNOWN WORKSPACE NUMBER"
    ;;
esac

# Create cromwell config file
if [ -e /home/jupyter/.cromwell ]; then
  rm -rf /home/jupyter/.cromwell
fi
mkdir /home/jupyter/.cromwell
wb cromwell generate-config \
  --google-bucket-name=$WORKSPACE_BUCKET \
  --dir=/home/jupyter/.cromwell
cat << EOF > /home/jupyter/.cromwell/cromwell.override.conf
include "cromwell.conf"

backend.providers.GCPBATCH.config {
  concurrent-job-limit = 100
}

call-caching {
  enabled = true
}
EOF

# Create cromshell config file
cat << EOF > /home/jupyter/.cromshell/cromshell_config.json
{
  "cromwell_server": "http://localhost:8000",
  "requests_timeout": 5
}
EOF

# Launch cromwell in server mode
java \
  -Xms8G \
  -Xmx16G \
  -Dconfig.file=/home/jupyter/.cromwell/cromwell.override.conf \
  -jar $CROMWELL_JAR \
  server

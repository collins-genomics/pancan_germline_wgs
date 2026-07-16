#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Cromwell server launch script for AoU RW v2.0 (Verily Pre)

# Note that this must be run in a separate terminal

# Create cromwell config file
mkdir /home/jupyter/.cromwell
wb cromwell generate-config \
  --google-bucket-name=$WORKSPACE_BUCKET \
  --dir=/home/jupyter/.cromwell

# Create cromshell config file
cat << EOF > /home/jupyter/.cromshell/cromshell_config.json
{
  "cromwell_server": "http://localhost:8000",
  "requests_timeout": 5
}
EOF

# Launch cromwell in server mode
java \
  -Dconfig.file=/home/jupyter/.cromwell/cromwell.conf \
  -jar $CROMWELL_JAR \
  server

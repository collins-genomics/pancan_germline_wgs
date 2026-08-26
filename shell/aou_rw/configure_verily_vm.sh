#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# VM configuration startup script for AoU RW v2.0 (Verily Pre)

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://rw-migration-aou-rw-84a0039b

# Confirm $MAIN_WORKSPACE_BUCKET is registered in this workspace as a 
# referenced resource, and re-register if not
if [ $( wb resource list | fgrep main_g2c_bucket | wc -l ) -lt 1 ]; then
  wb resource add-ref gcs-bucket \
    --id main_g2c_bucket \
    --bucket-name rw-migration-aou-rw-84a0039b
fi

# Copy necessary code to local disk
gcloud storage rsync -r $MAIN_WORKSPACE_BUCKET/code ./code
find code/ -name "*.py" | xargs -I {} chmod a+x {}
find code/ -name "*.R" | xargs -I {} chmod a+x {}
find code/ -name "*.sh" | xargs -I {} chmod a+x {}

# Upgrade to Cromwell v92
if ! [ -e ~/bin/cromwell-92.jar ]; then
  mkdir ~/bin
  wget -P ~/bin/ https://github.com/broadinstitute/cromwell/releases/download/92/cromwell-92.jar
fi
export CROMWELL_JAR=~/bin/cromwell-92.jar

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Install necessary packages
. code/refs/install_packages.sh python R
if [ "$( echo $CONDA_DEFAULT_ENV )" != "g2c" ]; then
  source activate g2c
fi

# Create dependencies .zip for generic G2C workflow submissions
cd code/wdl/pancan_germline_wgs && \
zip -r g2c.dependencies.zip . && \
mv g2c.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for QC workflow submissions
cd code/wdl/pancan_germline_wgs/vcf-qc && \
zip qc.dependencies.zip *.wdl && \
mv qc.dependencies.zip ~/ && \
cd ~

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

# Set up expected directory structure for local cromshell execution
for dir in cromshell cromshell/inputs cromshell/submissions cromshell/progress \
           cromshell/job_ids; do
  if ! [ -e ~/$dir ]; then
    mkdir ~/$dir
  fi
done

# Download workspace-specific contig lists
gcloud storage rsync -r \
  gs://dfci-g2c-refs/hg38/contig_lists \
  ./contig_lists

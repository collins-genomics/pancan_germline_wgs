#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Extract desired subsets of the AoU CDR v8 srWGS short variant callset for G2C benchmarking

# Note that this code is designed to be run on a spark dataproc cluster in the AoU Researcher Workbench

# Note 2: this code was defined for the AoU RW v1.0, and it may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0

# Recommended dataproc cluster specs:
# Main node: 8CPU, 52GB RAM, 150GB disk
# 40 workers @ 4CPU, 16GB RAM, 150GB disk
# 200 preemptible workers @ 4CPU, 16GB RAM, 150GB disk

# Copy Hail curation script to local
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/code/scripts/extract_aou_srwgs_short_variants.py \
  ./
chmod a+x extract_aou_srwgs_short_variants.py

# Iterate over all 24 primary contigs and extract their short variants in serial
for k in $( seq 1 22 ) X Y; do
  contig="chr$k"
  ./extract_aou_srwgs_short_variants.py $contig
done

#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Prioritize relatives of G2C cancer cases from AoU to be included in phase 1 


# Note: this code was defined for the AoU RW v1.0, and it may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in data data/cram_paths cromshell cromshell/inputs cromshell/job_ids \
           cromshell/progress; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.R" | xargs -I {} chmod a+x {}

# Source .bashrc
. ~/code/refs/dotfiles/aou.rw.bashrc

# Localize & format sample ID lists and manifests
. code/refs/setup_sample_info.sh

# Download AoU precomputed kinship metrics
gsutil -m -u $GPROJECT cp \
  gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness.tsv \
  ./

# Pool list of all G2C cases
find sample_lists/ -name "*samples.list" | fgrep -v control \
| xargs -I {} cat {} | sort -V | uniq > all_g2c_cases.samples.list

# Pool list of all G2C samples processed to date
find sample_lists/ -name "*samples.list" | fgrep -v eligible \
| xargs -I {} cat {} | sort -V | uniq > all_g2c.samples.list

# Copy Noah's list of all UFC cases
gsutil -m cp \
  gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ufc_cohort.list \
  ./

# Select relatives
~/code/scripts/select_aou_relatives.R \
  --kinship relatedness.tsv \
  --case-ids all_g2c_cases.samples.list \
  --ufc-case-ids ufc_cohort.list \
  --processed-ids all_g2c.samples.list \
  --outfile sample_lists/relatives.samples.list

# Copy list of selected relatives to main storage bucket
gsutil -m cp \
  sample_lists/relatives.samples.list \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/sample_lists/


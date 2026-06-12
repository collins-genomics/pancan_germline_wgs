#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-HC and GATK-SV (module 01) on all AoU samples selected for inclusion in G2C pilot

# Note that this code is designed to be run inside the AoU Researcher Workbench
# See aou_bash_utils.sh for custom function definitions used below


# Note2 : this code was defined for the AoU RW v1.0, and it may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export BOTO_CONFIG=~/code/refs/dotfiles/dfci-g2c.aou-rw.boto.cfg
export BOTO_PATH=~/code/refs/dotfiles/
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in data data/cram_paths cromshell cromshell/inputs cromshell/job_ids \
           cromshell/progress; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}

# Install necessary packages
. code/refs/install_packages.sh python

# Source .bashrc and bash utility functions
. ~/code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh
. code/refs/aou_bash_utils.sh

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Localize & format sample ID lists and manifests
. code/refs/setup_sample_info.sh


###########
# GATK-HC #
###########

# Main loop for the submission and tracking of GATK-HC workflows
# Before launching, this code will check the status of each sample to determine eligibility
# If already successful, it will stage final outputs and clean up temporary files
# It will also update counts in the main sample progress tracker table
while read cancer; do
  check_status gatk-hc $cancer
  update_status_table gatk-hc
  cleanup_garbage
  submit_workflows gatk-hc $cancer
done < cancers.list

# Main loop for the submission and tracking of GATK-HC gVCF postprocessing
# Includes both reblocking and reheadering
# Before launching, this code will check the status of each sample to determine eligibility
# If already successful, it will stage final outputs and clean up temporary files
# It will also update counts in the main sample progress tracker table
while read cancer; do
  check_status gvcf-pp $cancer
  update_status_table gvcf-pp
  cleanup_garbage
  submit_workflows gvcf-pp $cancer
done < cancers.list


###########
# GATK-SV #
###########

# Main loop for the submission and tracking of GATK-SV workflows
# Before launching, this code will check the status of each sample to determine eligibility
# If already successful, it will stage final outputs and clean up temporary files
# It will also update counts in the main sample progress tracker table
while read cancer; do
  check_status gatk-sv $cancer
  update_status_table gatk-sv
  cleanup_garbage
  submit_workflows gatk-sv $cancer
done < cancers.list


################
# READ METRICS #
################

# Main loop for the submission and tracking of read metric collection workflows
# Before launching, this code will check the status of each sample to determine eligibility
# It will also update counts in the main sample progress tracker table
while read cancer; do
  check_status read-metrics $cancer
  update_status_table read-metrics
  cleanup_garbage
  submit_workflows read-metrics $cancer
done < cancers.list


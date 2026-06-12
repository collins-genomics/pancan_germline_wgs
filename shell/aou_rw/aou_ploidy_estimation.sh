#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-SV ploidy estimation & dosage bias scoring (module 02) on all AoU samples selected for G2C phase 1

# Note that this code is designed to be run inside the AoU Researcher Workbench
# See general_bash_utils.sh for custom function definitions used below


# Note 2: these functions were defined for the AoU RW v1.0
# They may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in cromshell cromshell/inputs cromshell/job_ids cromshell/progress \
           data data/batch_lists data/aou_ploidy_outputs; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
for suffix in py R; do
  find code/ -name "*.${suffix}" | xargs -I {} chmod a+x {}
done

# Install necessary packages
. code/refs/install_packages.sh python

# Source .bashrc and bash utility functions
. ~/code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Install necessary packages
. code/refs/install_packages.sh R

# Prep GATK-SV dependencies
cd code/wdl/gatk-sv && \
zip gatksv.dependencies.zip *.wdl && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Localize list of all G2C AoU samples
gsutil -m -u $GPROJECT cp \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/other/all_g2c.phase1.samples.list \
  data/


#########################
# RUN PLOIDY ESTIMATION #
#########################

# Randomly divide AoU samples into batches of ~500
code/scripts/evenSplitter.R \
  -L 500 \
  --shuffle \
  data/all_g2c.phase1.samples.list \
  data/batch_lists/AoU.ploidy_estimation.batch_

# Make list of all batches
find data/batch_lists/ -name "AoU.ploidy_estimation.batch_*" \
| xargs -I {} basename {} | sort -V \
> data/ploidy_batches.list

# Submit one workflow per batch while also checking the status of prior submissions
touch cromshell/progress/aou.ploidy_estimation.batch_status.tsv
while read bid; do

  # Check status of most recent submission, if any
  if [ -e cromshell/job_ids/$bid.ploidy_estimation.job_ids.list ]; then
    last_uuid=$( tail -n1 cromshell/job_ids/$bid.ploidy_estimation.job_ids.list )
    status=$( cromshell --no_turtle -t 120 -mc status $last_uuid \
              | jq .status | tr -d '"' | tr '[A-Z]' '[a-z]' | sed 's/[ ]+/_/g' )
  else
    status="not_started"
  fi
  echo -e "Status of batch $bid: $status"

  # Update status tracker
  fgrep -wv $bid cromshell/progress/aou.ploidy_estimation.batch_status.tsv \
  | cat - <( echo -e "$bid\t$status" ) | sort -Vk1,1 \
  >> cromshell/progress/aou.ploidy_estimation.batch_status.tsv2
  mv cromshell/progress/aou.ploidy_estimation.batch_status.tsv2 \
     cromshell/progress/aou.ploidy_estimation.batch_status.tsv

  # Only submit if status is not_started or failed
  if [ $status == "not_started" ] || [ $status == "failed" ] || [ $status == "aborted" ]; then

    echo -e "Submitting ploidy estimation workflow for batch $bid"

    # Format input .json
    export BATCH=$bid
    export SAMPLES=$( awk -v ORS=", " '{ print "\""$1"\"" }' data/batch_lists/$bid | sed 's/, $//' )
    export COUNTS=$( awk -v ORS=", " -v uri_base="$MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/aou/gatk-sv/coverage/" \
                         '{ print "\""uri_base$1".counts.tsv.gz\"" }' data/batch_lists/$bid | sed 's/, $//' )
    ~/code/scripts/envsubst.py \
      -i code/refs/json/aou.gatk_sv.02_evidence_qc.inputs.template.json \
    | sed 's/\t//g' | paste -s -d\ \
    > cromshell/inputs/$bid.ploidy_estimation.inputs.json

    # Submit job and add job ID to list of jobs for this workflow
    cmd="cromshell --no_turtle -t 120 -mc submit"
    cmd="$cmd --options-json code/refs/json/aou.cromwell_options.default.json"
    cmd="$cmd --dependencies-zip gatksv.dependencies.zip"
    cmd="$cmd code/wdl/gatk-sv/EvidenceQC.wdl"
    cmd="$cmd cromshell/inputs/$bid.ploidy_estimation.inputs.json"
    eval $cmd | jq .id | tr -d '"' \
    >> cromshell/job_ids/$bid.ploidy_estimation.job_ids.list
  fi
done < data/ploidy_batches.list


#################
# STAGE OUTPUTS #
#################

# Write list of batches to stage
# (Note that this needs to be done outside the while loop so input does not change)
awk '{ if ($2=="succeeded") print $1 }' \
  cromshell/progress/aou.ploidy_estimation.batch_status.tsv \
> cromshell/progress/aou.ploidy_estimation.batches_to_stage.list

# Stage each batch in serial
while read bid; do

  # Make directory
  if ! [ -e data/aou_ploidy_outputs/$bid ]; then
    mkdir data/aou_ploidy_outputs/$bid
  fi
    
  # For convenience, write list of samples from batch to batch directory
  sort -V data/batch_lists/$bid | awk -v OFS="\t" '{ print $1, "aou" }' \
  | cat <( echo -e "#sample\tcohort" ) - \
  > data/aou_ploidy_outputs/$bid/$bid.sample_info.tsv

  # Process each output in serial
  while read descrip uri; do

    # Process only the outputs we want to save
    case $descrip in

      # Download median coverage, WGD plot, and WGD .tsv for permanent local storage
      EvidenceQC.WGD_scores|EvidenceQC.WGD_dist|EvidenceQC.bincov_median)
        gsutil -m cp $uri data/aou_ploidy_outputs/$bid/
        ;;

      # Download & clean QC table
      EvidenceQC.qc_table)
        gsutil -m cp $uri data/aou_ploidy_outputs/$bid/$bid.evidence_qc_table.tsv
        code/scripts/clean_module02_qc_table.R \
          data/aou_ploidy_outputs/$bid/$bid.evidence_qc_table.tsv
        ;;

      # Download & process ploidy plots tarball (don't need to save everything)
      EvidenceQC.ploidy_plots)
        gsutil -m cp $uri ./
        tar -xzvf ./$( basename $uri ) --directory ./
        find ./ploidy_est/*png | xargs -I {} mv {} data/aou_ploidy_outputs/$bid/
        mv ./ploidy_est/sample_sex_assignments.txt.gz data/aou_ploidy_outputs/$bid/
        rm -rf ./ploidy_est ./$( basename $uri )
        ;;
    esac
    
  done < <( cromshell -t 120 --no_turtle -mc list-outputs \
              $( tail -n1 cromshell/job_ids/$bid.ploidy_estimation.job_ids.list ) \
            | sed 's/: /\t/g' )
  
done < cromshell/progress/aou.ploidy_estimation.batches_to_stage.list

# Once all batches have been processed, combine the QC tables across all batches
head -n1 $( find data/aou_ploidy_outputs/ -name "*.evidence_qc_table.tsv" | head -n1 ) \
| sed 's/^#ID/#cohort\tsample/' \
> data/aou_ploidy_outputs/dfci-g2c-ploidy-estimation.aou.merged_qc_table.tsv
find data/aou_ploidy_outputs/ -name "*.evidence_qc_table.tsv" \
| xargs -I {} cat {} | grep -ve '^#' | awk -v OFS="\t" '{ print "aou", $0 }' | sort -Vk1,1 -k2,2V \
>> data/aou_ploidy_outputs/dfci-g2c-ploidy-estimation.aou.merged_qc_table.tsv
gzip -f data/aou_ploidy_outputs/dfci-g2c-ploidy-estimation.aou.merged_qc_table.tsv

# Copy processed data to G2C staging bucket
gsutil -m cp -r \
  data/aou_ploidy_outputs \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/


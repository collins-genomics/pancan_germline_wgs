#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper bash functions for processing samples from the All of Us cohort 
# on the AoU Researcher Workbench


# Note: makes use of helper functions in general_bash_utils.sh, which must be
# sourced before this file is used


# Note 2: these functions were defined for the AoU RW v1.0
# They may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0


# Submit workflows for one cancer type and a single workflow 
# (gatk-hc, gatk-sv, gcnv-pp, or read-metrics)
submit_workflows() {
  # Check inputs
  if [ $# -ne 2 ]; then
    echo "Must specify (gatk-hc|gatk-sv|gvcf-pp|read-metrics) and cancer type as first two positional arguments"
    return
  else
    wflow=$1
    cancer_sub=$2
  fi

  # Set workflow-specific parameters
  wflow_lower=$( echo $wflow | sed 's/-/_/g' )
  wflow_nospace=$( echo $wflow_lower | sed 's/_//g' )
  status_tsv=~/cromshell/progress/$cancer_sub.$wflow_lower.sample_progress.tsv
  case $wflow in
    "gatk-hc")
      wdl=~/code/wdl/gatk-hc/haplotypecaller-gvcf-gatk4.wdl
      inputs_json_prefix=gatk_hc
      gate_width=60
      gate_timeout=60m
      workflow_name="GATK-HC"
      sid_cram_list=~/data/cram_paths/$cancer_sub.cram_paths.tsv
      ;;
    "gatk-sv")
      wdl=~/code/wdl/gatk-sv/GatherSampleEvidence.wdl
      inputs_json_prefix=gatk_sv_module_01
      gate_width=25
      gate_timeout=45m
      workflow_name="GATK-SV"
      sid_cram_list=~/data/cram_paths/$cancer_sub.cram_paths.tsv
      cd code/wdl/gatk-sv && \
      zip gatksv.dependencies.zip *.wdl && \
      mv gatksv.dependencies.zip ~/ && \
      cd ~
      ;;
    "gvcf-pp")
      wdl=~/code/wdl/pancan_germline_wgs/PostprocessGvcf.wdl
      inputs_json_prefix=gvcf_pp
      gate_width=150
      gate_timeout=30m
      workflow_name="gVCF postprocessing"
      awk '{ if ($2=="staged") print $1 }' \
        cromshell/progress/$cancer_sub.gatk_hc.sample_progress.tsv \
      | fgrep -wf - data/cram_paths/$cancer_sub.cram_paths.tsv \
      > $cancer_sub.$wflow.sids_to_submit.list
      sid_cram_list=$cancer_sub.$wflow.sids_to_submit.list
      ;;
    "read-metrics")
      wdl=~/code/wdl/pancan_germline_wgs/CalcReadPairProperties.wdl
      inputs_json_prefix=read_metrics
      gate_width=150
      gate_timeout=30m
      workflow_name="Read metric collection"
      sid_cram_list=~/data/cram_paths/$cancer_sub.cram_paths.tsv
      ;;
  esac
  if ! [ -e $status_tsv ]; then
    echo -e "Status tracker $status_tsv not found; skipping $workflow_name submissions for $cancer_sub patients"
    return 0
  fi

  # Only submit tasks for eligible samples (depending on last known status)
  k=0; j=0; s=0; g=0
  n=$( cat ~/data/cram_paths/$cancer_sub.cram_paths.tsv | wc -l )
  while read sid CRAM CRAI; do
    ((k++)); ((j++))
    # Check if sample has not been launched or has failed/aborted
    status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' $status_tsv )
    if [ -z $status ] || [ $status == "not_started" ] || \
       [ $status == "failed" ] || [ $status == "aborted" ] || \
       [ $status == "unknown" ]; then

      # Format sample-specific input .json
      ((s++)); ((g++))
      export sid=$sid
      export SID=$sid
      export CRAM=$CRAM
      export CRAI=$CRAI
      ~/code/scripts/envsubst.py \
        -i code/refs/json/aou.$inputs_json_prefix.inputs.template.json \
      | sed 's/\t//g' | paste -s -d\ \
      > cromshell/inputs/$sid.$wflow_nospace.inputs.json

      # Submit job and add job ID to list of jobs for this sample
      cmd="cromshell --no_turtle -t 120 -mc submit"
      cmd="$cmd --options-json code/refs/json/aou.cromwell_options.default.json"
      if [ $wflow == "gatk-sv" ]; then
        cmd="$cmd --dependencies-zip gatksv.dependencies.zip"
      fi
      cmd="$cmd $wdl cromshell/inputs/$sid.$wflow_nospace.inputs.json"
      eval $cmd | jq .id | tr -d '"' \
      >> cromshell/job_ids/$sid.$wflow_lower.job_ids.list
    fi

    # Only report on progress every 50 samples to limit terminal verbosity
    # This is one of the possible causes of large data egress warnings
    if [ $j -eq 50 ] || [ $k -eq $n ]; then
      echo -e "$k of $n $cancer_sub samples evaluated; $s $workflow_name jobs submitted"
      j=0
    fi

    # Only submit $gate_width workflows every $gate_timeout at most
    # This is in the hopes of limiting data egress surges to not trigger security warnings
    if [ $g -eq $gate_width ]; then
      echo -e "$g $workflow_name workflows submitted recently; waiting $gate_timeout before resuming submissions..."
      g=0
      sleep $gate_timeout
    fi
  done < $sid_cram_list
  if [ $wflow == "gvcf-pp" ]; then
    rm $cancer_sub.$wflow.sids_to_submit.list
  fi
}


# Check sample status for one cancer type and a single workflow 
# (gatk-hc, gatk-sv, gcnv-pp, or read-metrics)
check_status() {
  if [ $# -ne 2 ]; then
    echo "Must specify (gatk-hc|gatk-sv|gvcf-pp|read-metrics) and cancer type as first two positional arguments"
    return
  else
    wflow=$1
    cancer_sub=$2
  fi
  status_tsv=~/cromshell/progress/$cancer_sub.$( echo $wflow | sed 's/-/_/g' ).sample_progress.tsv
  if [ -e $status_tsv ]; then
    n=$( cat ~/sample_lists/$cancer_sub.samples.list | wc -l )
    k=0
    j=0
    while read sid cram crai; do
      ((k++))
      ((j++))
      # Only print to stdout every 50 samples to limit terminal verbosity
      # This is one of the possible causes of large data egress warnings
      if [ $j -eq 50 ]; then
        echo -e "Checking $wflow for $sid; sample $k of $n $cancer_sub patients"
        j=0
      fi
      ~/code/scripts/check_aou_gatk_status.py \
        --sample-id "$sid" \
        --mode $wflow \
        --bucket $WORKSPACE_BUCKET \
        --staging-bucket $MAIN_WORKSPACE_BUCKET \
        --status-tsv $status_tsv \
        --update-status \
        --unsafe \
        --metrics-optional \
      2> /dev/null
    done < ~/data/cram_paths/$cancer_sub.cram_paths.tsv
  fi
}


# Update status table of sample progress
update_status_table() {
  if [ $# -ne 1 ]; then
    echo "Must specify (gatk-hc|gatk-sv|gvcf-pp|read-metrics) as first positional argument"
    return
  else
    wflow=$1
  fi
  while read cancer_sub; do
    cancer_progress_tsv=~/cromshell/progress/$cancer_sub.$( echo $wflow | sed 's/-/_/g' ).sample_progress.tsv
    if ! [ -e $cancer_progress_tsv ]; then
      continue
    fi
    case $cancer_sub in
      *control*)
        name_for_table="controls"
        ;;
      pediatric)
        name_for_table="other"
        ;;
      *)
        name_for_table=$cancer_sub
        ;;
    esac
    awk -v cancer=$name_for_table -v wflow=$wflow \
      '{ if ($1!=cancer || $2!=wflow) print $0 }' \
      ~/cromshell/progress/gatk.sample_progress.summary.tsv \
    > ~/cromshell/progress/gatk.sample_progress.summary.tsv2
    cut -f2 $cancer_progress_tsv \
    | sort | uniq -c | sort -nrk1,1 \
    | awk -v cancer=$name_for_table -v wflow=$wflow -v OFS="\t" \
      '{ print cancer, wflow, $2, $1 }' \
    >> ~/cromshell/progress/gatk.sample_progress.summary.tsv2
    mv \
      ~/cromshell/progress/gatk.sample_progress.summary.tsv2 \
      ~/cromshell/progress/gatk.sample_progress.summary.tsv
  done < ~/cancers.list
  sort -Vk1,3 ~/cromshell/progress/gatk.sample_progress.summary.tsv
}


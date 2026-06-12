#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper bash functions for running the GATK-SV cohort mode pipeline 
# on the All of Us Researcher Workbench


# Note: makes use of helper functions in general_bash_utils.sh, which must be
# sourced before this file is used


# Note 2: these functions were defined for the AoU RW v1.0
# They may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0


# Check the status of a GATK-SV module for a single batch
# Also updates running tracker
# Also manages submission or staging of outputs as needed
check_batch_module() {
   # Check inputs
  if [ $# -ne 2 ]; then
    echo "Must specify batch ID and [03-08,10,14A] as first two positional arguments"
    return
  else
    export bid=$1
    export module_idx=$2
  fi

  # Check to confirm tracker exists; if it doesn't, create it
  tracker="cromshell/progress/gatksv.batch_modules.progress.tsv"
  if ! [ -e $tracker ]; then touch $tracker; fi

  # Update status of most recent submission for batch
  ~/code/scripts/check_gatksv_batch_module_status.py \
    --batch-id "$bid" \
    --module-index "$module_idx" \
    --staging-bucket "$MAIN_WORKSPACE_BUCKET" \
    --status-tsv "$tracker" \
    --update-status

  # Check number of prior submissions (this can impact gate window)
  jid_list=~/cromshell/job_ids/$bid.$sub_name.job_ids.list
  if [ -e $jid_list ]; then
    n_prev_subs=$( cat $jid_list | wc -l )
  else
    n_prev_subs=0
  fi

  # Get module-specific parameters
  case $module_idx in
    03)
      sub_name="03-TrainGCNV"
      gate=1
      max_resub=3
      ;;
    04)
      sub_name="04-GatherBatchEvidence"
      gate=30
      max_resub=15
      ;;
    05)
      sub_name="05-ClusterBatch"
      gate=2
      max_resub=2
      ;;
    05B)
      sub_name="05B-ExcludeClusteredOutliers"
      gate=0
      max_resub=1
      ;;
    05C)
      sub_name="05C-ReclusterBatch"
      gate=2
      max_resub=2
      ;;
    06)
      sub_name="06-GenerateBatchMetrics"
      gate=30
      max_resub=5
      ;;
    07)
      sub_name="07-FilterBatchSites"
      gate=20
      max_resub=2
      ;;
    08)
      sub_name="08-FilterBatchSamples"
      gate=1
      max_resub=2
      ;;
    10)
      sub_name="10-GenotypeBatch"
      if [ $n_prev_subs -gt 0 ]; then
        gate=5
      else
        gate=60
      fi
      max_resub=3
      ;;
    14A)
      sub_name="14A-FilterCoverageSamples"
      gate=1
      max_resub=3
      ;;
  esac

  # Process batch based on reported status
  status=$( awk -v bid="$bid" -v midx="$module_idx" \
              '{ if ($1==bid && $2==midx) print $3 }' \
              $tracker )
  if [ -z $status ]; then
    echo -e "Failed getting status of batch '$bid' for module '$module_idx'. Exiting."
    return 1
  fi
  case "$status" in
    not_started|failed|doomed|aborted)
      if [ $n_prev_subs -le $max_resub ]; then
        echo -e "Submitting GATK-SV module $module_idx for batch $bid"
        submit_batch_module "$bid" "$module_idx"
        echo -e "Waiting $gate minutes before continuing to avoid choking Cromwell server..."
        sleep ${gate}m
      else
        echo -e "Batch $bid has hit the maximum of $max_resub submissions for GATK-SV module $module_idx. Not submitting a new workflow for this batch."
      fi
      ;;
  esac
}


# Loop to manage submissions & progress tracking for a GATK-SV module for all batches
module_submission_routine_all_batches() {
   # Check input
  if [ $# -ne 1 ]; then
    echo "Must specify [03-08,10] as only positional argument"
    return
  else
    export module_idx=$1
  fi

  WN=$( get_workspace_number )
  tracker=cromshell/progress/gatksv.batch_modules.progress.tsv

  # Customize gate duration for certain short-running workflows
  outer_gate=60
  case $module_idx in
    05|05C)
      outer_gate=30
      ;;
    05B|08)
      outer_gate=10
      ;;
    07)
      outer_gate=20
      ;;
  esac

  # Batches will always be parallelized by workspace except for 14A
  if [ $module_idx == "14A" ]; then
    batches_list="batch_info/dfci-g2c.gatk-sv.batches.list"
  else
    batches_list="batch_info/dfci-g2c.gatk-sv.batches.$WN.list"
  fi

  _count_remaining() {
    if [ -e $tracker ]; then
      incomplete=$( awk -v midx=$1 '{ if ($2==midx && $3!="staged") print }' \
                      $tracker | wc -l )
      missing=$( awk -v midx=$1 '{ if ($2==midx) print $1 }' $tracker \
                 | fgrep -wvf - $batches_list | wc -l )
    else
      incomplete=0
      missing=$( cat $batches_list | wc -l )
    fi
    echo $(( $incomplete + $missing ))
  }

  while [ $( _count_remaining $module_idx ) -gt 0 ]; do
    while read bid; do
      echo -e "Checking status of $bid for GATK-SV module $module_idx"
      check_batch_module $bid $module_idx
      cleanup_garbage
    done < $batches_list
    echo -e "Finished checking status of all batches for GATK-SV module $module_idx"
    echo -e "Status of all batches for all GATK-SV modules:"
    report_gatksv_status $module_idx
    if [ $( _count_remaining $module_idx ) -eq 0 ]; then
      break
    fi
    echo -e "Waiting $outer_gate minutes before checking progress again..."
    sleep ${outer_gate}m
  done

  echo -e "All batches finished for GATK-SV module $module_idx. Ending monitor routine."
}


# Helper function to report progress for a single module for all batches
report_gatksv_status() {
  # Optional first argument: module index
  tracker="cromshell/progress/gatksv.batch_modules.progress.tsv"
  if [ -e $tracker ]; then
    if [ $# -eq 0 ]; then
      cat $tracker
    else
      awk -v query=$1 -v FS="\t" -v OFS="\t" '{ if ($2==query) print }' $tracker
    fi
  fi
}


# Submit a single GATK-SV module for a single batch
submit_batch_module() {
   # Check inputs
  if [ $# -ne 2 ]; then
    echo "Must specify batch ID and [03-08,10,14A] as first two positional arguments"
    return
  else
    export BATCH=$1
    export module_idx=$2
  fi

  # Set constants
  g2c_mfst="dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz"
  batch_sid_list="batch_info/sample_lists/$BATCH.samples.list"

  # For modules after 03, check to ensure that the previous module has been staged
  case $module_idx in
    05B)
      prev_idx="04"
      ;;
    05C)
      prev_idx="05B"
      ;;
    06)
      prev_idx="05C"
      ;;
    10)
      prev_idx="08"
      ;;
    14A)
      prev_idx="04"
      ;;
    *)
      prev_idx="0$(( 10#$module_idx - 1 ))"    
      ;;
  esac
  prev_module_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/$prev_idx/$BATCH/$BATCH.gatksv_module_$prev_idx.outputs.json
  if [ $( gsutil ls $prev_module_outputs_json | wc -l ) -lt 1 ]; then
    echo "Module $module_idx requires the outputs from module $prev_idx to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
    return 126
  fi

  # Some modules require multiple prior outputs; these extra outputs are checked here
  case $module_idx in
    06)
      # Also requires 04 output
      module_04_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/04/$BATCH/$BATCH.gatksv_module_04.outputs.json
      if [ $( gsutil ls $module_04_outputs_json | wc -l ) -lt 1 ]; then
        echo "Module $module_idx requires the outputs from module 04 to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
        return 126
      fi
      ;;
    07)
      # Also requires 05C output
      module_05C_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/05C/$BATCH/$BATCH.gatksv_module_05C.outputs.json
      if [ $( gsutil ls $module_05C_outputs_json | wc -l ) -lt 1 ]; then
        echo "Module $module_idx requires the outputs from module 05C to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
        return 126
      fi
      ;;
    10)
      # Also requires 04 output
      module_04_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/04/$BATCH/$BATCH.gatksv_module_04.outputs.json
      if [ $( gsutil ls $module_04_outputs_json | wc -l ) -lt 1 ]; then
        echo "Module $module_idx requires the outputs from module 04 to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
        return 126
      fi
      # Also requires 07 output
      module_07_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/07/$BATCH/$BATCH.gatksv_module_07.outputs.json
      if [ $( gsutil ls $module_07_outputs_json | wc -l ) -lt 1 ]; then
        echo "Module $module_idx requires the outputs from module 07 to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
        return 126
      fi
      ;;
  esac

  # Make batch staging directory, if necessary
  for dir in staging "staging/$BATCH"; do
    if ! [ -e $dir ]; then
      mkdir ~/$dir
    fi
  done
  case $module_idx in
    03)
      module_name="TrainGCNV"
      ;;
    04)
      module_name="GatherBatchEvidence"
      ;;
    05)
      module_name="ClusterBatch"
      ;;
    05B)
      module_name="ExcludeClusteredOutliers"
      ;;
    05C)
      module_name="ReclusterBatch"
      ;;
    06)
      module_name="GenerateBatchMetrics"
      ;;
    07)
      module_name="FilterBatchSites"
      ;;
    08)
      module_name="FilterBatchSamples"
      ;;
    10)
      module_name="GenotypeBatch"
      ;;
    14A)
      module_name="FilterCoverageSamples"
      ;;
    *)
      echo "Module number $module_idx not recognized by submit_batch_module. Exiting."
      return 2
      ;;
  esac
  sub_name="${module_idx}-$module_name"
  sub_dir="staging/$BATCH/$sub_name"
  if ! [ -e $sub_dir ]; then mkdir $sub_dir; fi

  # Get map of G2C ID, original ID, and original cohort for this batch
  zcat $g2c_mfst | fgrep -wf $batch_sid_list | cut -f1-3 \
  > staging/$BATCH/$batch.sample_info.tsv
  cut -f1 staging/$BATCH/$batch.sample_info.tsv \
  > staging/$BATCH/$batch.samples.list

  # Set workflow-specific parameters
  json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.$sub_name.inputs.template.json
  case $module_idx in

    #############
    # MODULE 03 #
    #############
    03)
      wdl="code/wdl/gatk-sv/TrainGCNV.wdl"
      while read sid oid cohort; do
        if [ "$cohort" == "aou" ]; then
          gs_base="$MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs"
        else
          gs_base="gs://dfci-g2c-inputs"
        fi
        echo "$gs_base/$cohort/gatk-sv/coverage/$oid.counts.tsv.gz"
      done < staging/$BATCH/$batch.sample_info.tsv \
      > $sub_dir/$BATCH.cov.list
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "GatherBatchEvidence.counts" : $( collapse_txt $sub_dir/$BATCH.cov.list ),
    "GatherBatchEvidence.samples" : $( collapse_txt staging/$BATCH/$batch.samples.list )
}
EOF
      ;;

    #############
    # MODULE 04 #
    #############
    04)
      wdl="code/wdl/gatk-sv/GatherBatchEvidence.wdl"

      # Build necessary input arrays
      for suf in cov.list pe.list sr.list sd.list manta.list melt.list wham.list; do
        lfile="$sub_dir/${BATCH}.$suf"
        if [ -e $lfile ]; then rm $lfile; fi
      done
      while read sid oid cohort; do
        if [ "$cohort" == "aou" ]; then
          gs_base="$MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs"
        else
          gs_base="gs://dfci-g2c-inputs"
        fi
        echo "$gs_base/$cohort/gatk-sv/coverage/$oid.counts.tsv.gz" >> $sub_dir/$BATCH.cov.list
        echo "$gs_base/$cohort/gatk-sv/pesr/$oid.pe.txt.gz" >> $sub_dir/$BATCH.pe.list
        echo "$gs_base/$cohort/gatk-sv/pesr/$oid.sr.txt.gz" >> $sub_dir/$BATCH.sr.list
        echo "$gs_base/$cohort/gatk-sv/pesr/$oid.sd.txt.gz" >> $sub_dir/$BATCH.sd.list
        echo "$gs_base/$cohort/manta/$oid.manta.vcf.gz" >> $sub_dir/$BATCH.manta.list
        echo "$gs_base/$cohort/melt/$oid.melt.vcf.gz" >> $sub_dir/$BATCH.melt.list
        echo "$gs_base/$cohort/wham/$oid.wham.vcf.gz" >> $sub_dir/$BATCH.wham.list
      done < staging/$BATCH/$batch.sample_info.tsv

      # Begin building .json template updates
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "GatherBatchEvidence.PE_files" : $( collapse_txt $sub_dir/$BATCH.pe.list ),
    "GatherBatchEvidence.SD_files" : $( collapse_txt $sub_dir/$BATCH.sd.list ),
    "GatherBatchEvidence.SR_files" : $( collapse_txt $sub_dir/$BATCH.sr.list ),
    "GatherBatchEvidence.contig_ploidy_model_tar" : $( gsutil cat $prev_module_outputs_json \
                                                       | jq .cohort_contig_ploidy_model_tar \
                                                       | awk '{ if ($1 ~ /gs/) print $1 }' ),
    "GatherBatchEvidence.counts" : $( collapse_txt $sub_dir/$BATCH.cov.list ),
    "GatherBatchEvidence.gcnv_model_tars" : $( gsutil cat $prev_module_outputs_json \
                                               | jq .cohort_gcnv_model_tars \
                                               | paste -s - -d\ | tr -s ' ' ),
    "GatherBatchEvidence.manta_vcfs" : $( collapse_txt $sub_dir/$BATCH.manta.list ),
    "GatherBatchEvidence.melt_vcfs" : $( collapse_txt $sub_dir/$BATCH.melt.list ),
EOF

      # Get smart about resubmissions: try to determine failure mode and 
      # increase memory for just that failing job
      if [ -e cromshell/job_ids/$BATCH.$sub_name.job_ids.list ]; then
        n_prev_subs=$( cat cromshell/job_ids/$BATCH.$sub_name.job_ids.list | wc -l )
      else
        n_prev_subs=0
      fi
      if [ $n_prev_subs -gt 0 ]; then
        last_sub_id=$( tail -n1 cromshell/job_ids/$BATCH.$sub_name.job_ids.list )
        check_cromwell_return_codes \
          $WORKSPACE_BUCKET/cromwell-execution/$module_name/$last_sub_id \
        > $sub_dir/$BATCH.$last_sub_id.fail_rcs.list
  
        # gCNV case mode
        gcnv_ram=$( jq '."GatherBatchEvidence.runtime_attr_case"' \
                      cromshell/inputs/$BATCH.$sub_name.inputs.json \
                    | fgrep mem_gb | awk '{ print $NF }' )
        gcnv_fails=$( fgrep gCNV $sub_dir/$BATCH.$last_sub_id.fail_rcs.list | wc -l )
        if [ $gcnv_fails -gt 0 ] || [ ! -z $gcnv_ram ]; then
          # Get previous cn.MOPS NormalR1 memory
          if [ -z $gcnv_ram ]; then
            gcnv_ram=48
          fi
          if [ $gcnv_fails -gt 0 ]; then
            # Increment prior memory by +4 GB up to max permitted of 64GB
            gcnv_ram=$(( $gcnv_ram + 4 ))
          fi
          if [ $gcnv_ram -gt 32 ]; then
            gcnv_ram=32
            echo -e "WARNING: module $module_idx for $BATCH is requesting more memory for GATK-gCNV but has already been attempted at max allowed value (32GB)"
          fi
          cat << EOF >> $sub_dir/$BATCH.$sub_name.updates.json
    "GatherBatchEvidence.runtime_attr_case" : { "mem_gb" : $gcnv_ram },
EOF
        fi

        # cn.MOPS
        cnmops_ram=$( jq '."GatherBatchEvidence.cnmops_sample3_runtime_attr"' \
                          cromshell/inputs/$BATCH.$sub_name.inputs.json \
                    | fgrep mem_gb | awk '{ print $NF }' )
        cnmops_fails=$( fgrep CNMOPS $sub_dir/$BATCH.$last_sub_id.fail_rcs.list | wc -l )
        if [ $cnmops_fails -gt 0 ] || [ ! -z $cnmops_ram ]; then
          # Get previous cn.MOPS NormalR1 memory
          if [ -z $cnmops_ram ]; then
            cnmops_ram=48
          fi
          if [ $cnmops_fails -gt 0 ]; then
            # Increment prior memory by +4 GB up to max permitted of 64GB
            cnmops_ram=$(( $cnmops_ram + 4 ))
          fi
          if [ $cnmops_ram -gt 64 ]; then
            cnmops_ram=64
            echo -e "WARNING: module $module_idx for $BATCH is requesting more memory for cn.MOPS but has already been attempted at max allowed value (64GB)"
          fi
          cat << EOF >> $sub_dir/$BATCH.$sub_name.updates.json
    "GatherBatchEvidence.cnmops_sample3_runtime_attr" : { "mem_gb" : $cnmops_ram , "boot_disk_gb" : 20, "disk_gb" : 100 },
EOF
        fi
      fi

      # Finish building .json template updates
      cat << EOF >> $sub_dir/$BATCH.$sub_name.updates.json
    "GatherBatchEvidence.samples" : $( collapse_txt staging/$BATCH/$batch.samples.list ),
    "GatherBatchEvidence.wham_vcfs" : $( collapse_txt $sub_dir/$BATCH.wham.list )
}
EOF
      ;;

    #############
    # MODULE 05 #
    #############
    05)
      wdl="code/wdl/gatk-sv/ClusterBatch.wdl"
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "ClusterBatch.del_bed" : $( gsutil cat $prev_module_outputs_json | jq .merged_dels ),
    "ClusterBatch.dup_bed" : $( gsutil cat $prev_module_outputs_json | jq .merged_dups ),
    "ClusterBatch.manta_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_manta_vcf_tar ),
    "ClusterBatch.melt_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_melt_vcf_tar ),
    "ClusterBatch.wham_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_wham_vcf_tar )
}
EOF
      ;;

    ##############
    # MODULE 05B #
    ##############
    05B)
      wdl="code/wdl/pancan_germline_wgs/ExcludeClusteredOutliers.wdl"
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "ExcludeClusteredOutliers.del_bed" : $( gsutil cat $prev_module_outputs_json | jq .merged_dels ),
    "ExcludeClusteredOutliers.dup_bed" : $( gsutil cat $prev_module_outputs_json | jq .merged_dups ),
    "ExcludeClusteredOutliers.manta_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_manta_vcf_tar ),
    "ExcludeClusteredOutliers.melt_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_melt_vcf_tar ),
    "ExcludeClusteredOutliers.wham_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_wham_vcf_tar )
}
EOF
      ;;

    ##############
    # MODULE 05C #
    ##############
    05C)
      wdl="code/wdl/gatk-sv/ClusterBatch.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.05-ClusterBatch.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "ClusterBatch.del_bed" : $( gsutil cat $prev_module_outputs_json | jq .del_bed_cleaned ),
    "ClusterBatch.dup_bed" : $( gsutil cat $prev_module_outputs_json | jq .dup_bed_cleaned ),
    "ClusterBatch.manta_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .manta_vcf_tar_cleaned ),
    "ClusterBatch.melt_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .melt_vcf_tar_cleaned ),
    "ClusterBatch.wham_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .wham_vcf_tar_cleaned )
}
EOF
      ;;

    #############
    # MODULE 06 #
    #############
    06)
      wdl="code/wdl/gatk-sv/GenerateBatchMetrics.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.06-GenerateBatchMetrics.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "GenerateBatchMetrics.baf_metrics" : $( gsutil cat $module_04_outputs_json | jq .merged_BAF ),
    "GenerateBatchMetrics.coveragefile" : $( gsutil cat $module_04_outputs_json | jq .merged_bincov ),
    "GenerateBatchMetrics.depth_vcf" : $( gsutil cat $prev_module_outputs_json | jq .clustered_depth_vcf ),
    "GenerateBatchMetrics.discfile" : $( gsutil cat $module_04_outputs_json | jq .merged_PE ),
    "GenerateBatchMetrics.manta_vcf" : $( gsutil cat $prev_module_outputs_json | jq .clustered_manta_vcf ),
    "GenerateBatchMetrics.medianfile" : $( gsutil cat $module_04_outputs_json | jq .median_cov ),
    "GenerateBatchMetrics.melt_vcf" : $( gsutil cat $prev_module_outputs_json | jq .clustered_melt_vcf ),
    "GenerateBatchMetrics.splitfile" : $( gsutil cat $module_04_outputs_json | jq .merged_SR ),
    "GenerateBatchMetrics.wham_vcf" : $( gsutil cat $prev_module_outputs_json | jq .clustered_wham_vcf )
}
EOF
      ;;

    #############
    # MODULE 07 #
    #############
    07)
      wdl="code/wdl/gatk-sv/FilterBatchSites.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.07-FilterBatchSites.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "FilterBatchSites.depth_vcf": $( gsutil cat $module_05C_outputs_json | jq .clustered_depth_vcf ),
    "FilterBatchSites.evidence_metrics": $( gsutil cat $prev_module_outputs_json | jq .metrics ),
    "FilterBatchSites.evidence_metrics_common": $( gsutil cat $prev_module_outputs_json | jq .metrics_common ),
    "FilterBatchSites.manta_vcf": $( gsutil cat $module_05C_outputs_json | jq .clustered_manta_vcf ),
    "FilterBatchSites.melt_vcf": $( gsutil cat $module_05C_outputs_json | jq .clustered_melt_vcf ),
    "FilterBatchSites.wham_vcf": $( gsutil cat $module_05C_outputs_json | jq .clustered_wham_vcf )
}
EOF
      ;;

    #############
    # MODULE 08 #
    #############
    08)
      wdl="code/wdl/gatk-sv/FilterBatchSamples.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.08-FilterBatchSamples.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "FilterBatchSamples.depth_vcf": $( gsutil cat $prev_module_outputs_json | jq .sites_filtered_depth_vcf ),
    "FilterBatchSamples.manta_vcf": $( gsutil cat $prev_module_outputs_json | jq .sites_filtered_manta_vcf ),
    "FilterBatchSamples.melt_vcf": $( gsutil cat $prev_module_outputs_json | jq .sites_filtered_melt_vcf ),
    "FilterBatchSamples.wham_vcf": $( gsutil cat $prev_module_outputs_json | jq .sites_filtered_wham_vcf )
}
EOF
      ;;

    #############
    # MODULE 10 #
    #############
    10)
      wdl="code/wdl/gatk-sv/GenotypeBatch.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.10-GenotypeBatch.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "GenotypeBatch.batch_depth_vcf": $( gsutil cat $prev_module_outputs_json | jq .outlier_filtered_depth_vcf ),
    "GenotypeBatch.batch_pesr_vcf": $( gsutil cat $prev_module_outputs_json | jq .outlier_filtered_pesr_vcf ),
    "GenotypeBatch.coveragefile": $( gsutil cat $module_04_outputs_json | jq .merged_bincov ),
    "GenotypeBatch.discfile": $( gsutil cat $module_04_outputs_json | jq .merged_PE ),
    "GenotypeBatch.medianfile": $( gsutil cat $module_04_outputs_json | jq .median_cov ),
    "GenotypeBatch.rf_cutoffs": $( gsutil cat $module_07_outputs_json | jq .cutoffs ),
    "GenotypeBatch.splitfile": $( gsutil cat $module_04_outputs_json | jq .merged_SR )
}
EOF
      ;;

    ##############
    # MODULE 14A #
    ##############
    14A)
      wdl="code/wdl/gatk-sv/FilterCoverageSamples.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.14A-FilterCoverageSamples.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "FilterCoverageSamples.median_cov" : $( gsutil cat $prev_module_outputs_json | jq .median_cov )
}
EOF
      ;;

    *)
      echo "Module number $module_idx not recognized by submit_batch_module. Exiting."
      return 2
      ;;
  
  esac

  # Format input .json
  ~/code/scripts/envsubst.py \
    -i $json_input_template \
  | ~/code/scripts/update_json.py \
    -i stdin \
    -u $sub_dir/$BATCH.$sub_name.updates.json \
    -o cromshell/inputs/$BATCH.$sub_name.inputs.json

  # Submit job and add job ID to list of jobs for this batch
  cmd="cromshell --no_turtle -t 120 -mc submit"
  cmd="$cmd --options-json code/refs/json/aou.cromwell_options.default.json"
  cmd="$cmd --dependencies-zip gatksv.dependencies.zip"
  cmd="$cmd $wdl cromshell/inputs/$BATCH.$sub_name.inputs.json"
  eval $cmd | jq .id | tr -d '"' \
  >> cromshell/job_ids/$BATCH.$sub_name.job_ids.list
}


# Submit a single GATK-SV module for the entire cohort
submit_cohort_module() {
   # Check inputs
  if [ $# -ne 1 ]; then
    echo "Must specify [09,11-19] as first two positional arguments"
    return
  else
    export module_idx=$1
  fi

  # Set constants
  batches_list="batch_info/dfci-g2c.gatk-sv.batches.list"

  # Make batch staging directory, if necessary
  if ! [ -e ~/staging ]; then mkdir ~/staging; fi
  case $module_idx in
    09)
      module_name="MergeBatchSites"
      ;;
    11)
      module_name="RegenotypeCNVs"
      ;;
    12)
      module_name="CombineBatches"
      max_attempts=5
      ;;
    13)
      module_name="ResolveComplexVariants"
      ;;
    14)
      module_name="GenotypeComplexVariants"
      ;;
    15)
      module_name="CleanVcf"
      ;;
    16)
      module_name="RefineComplexVariants"
      max_attempts=6
      ;;
    17)
      module_name="JoinRawCalls"
      ;;
    18)
      module_name="SVConcordance"
      max_attempts=3
      ;;
    19)
      module_name="FilterGenotypes"
      max_attempts=2
      ;;
    *)
      echo "Module number $module_idx not recognized by submit_cohort_module. Exiting."
      return 2
      ;;
  esac
  sub_name="${module_idx}-$module_name"
  sub_dir=~/staging/$sub_name
  if [ -e $sub_dir ]; then rm -rf $sub_dir; fi
  mkdir $sub_dir

  # Set workflow-specific parameters
  staging_prefix=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs
  case $module_idx in

    #############
    # MODULE 09 #
    #############
    09)
      wdl="code/wdl/gatk-sv/MergeBatchSites.wdl"

      # Get URIs per batch
      while read bid; do
        mod08_output_json=$staging_prefix/08/$bid/$bid.gatksv_module_08.outputs.json
        gsutil cat $mod08_output_json | jq .outlier_filtered_depth_vcf \
        | tr -d '"' >> $sub_dir/depth_vcfs.list
        gsutil cat $mod08_output_json | jq .outlier_filtered_pesr_vcf \
        | tr -d '"' >> $sub_dir/pesr_vcfs.list
      done < $batches_list

      # Prep input .json
      cat << EOF > cromshell/inputs/dfci-g2c.v1.$sub_name.inputs.json
{
    "MergeBatchSites.cohort": "dfci-g2c.v1",
    "MergeBatchSites.depth_vcfs": $( collapse_txt $sub_dir/depth_vcfs.list ),
    "MergeBatchSites.pesr_vcfs": $( collapse_txt $sub_dir/pesr_vcfs.list ),
    "MergeBatchSites.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-11-15-v1.0-488d7cb0",
    "MergeBatchSites.runtime_attr_merge_pesr" : { "disk_gb" : 25 }
}
EOF
      ;;

    #############
    # MODULE 11 #
    #############
    11)
      wdl="code/wdl/gatk-sv/RegenotypeCNVs.wdl"

      # Collect batch-specific URIs and other info
      while read bid; do
        
        mod04_output_json=$staging_prefix/04/$bid/$bid.gatksv_module_04.outputs.json
        gsutil cat $mod04_output_json | jq .merged_bincov | tr -d '"' \
        >> $sub_dir/bincovs.list
        gsutil cat $mod04_output_json | jq .merged_bincov_index | tr -d '"' \
        >> $sub_dir/bincov_idxs.list
        gsutil cat $mod04_output_json | jq .median_cov | tr -d '"' \
        >> $sub_dir/cov_medians.list

        mod08_output_json=$staging_prefix/08/$bid/$bid.gatksv_module_08.outputs.json
        gsutil cat $mod08_output_json | jq .outlier_filtered_depth_vcf \
        | tr -d '"' >> $sub_dir/depth_vcfs.list
        
        mod10_output_json=$staging_prefix/10/$bid/$bid.gatksv_module_10.outputs.json
        gsutil cat $mod10_output_json | jq .genotyped_depth_vcf | tr -d '"' \
        >> $sub_dir/depth_gt_vcfs.list
        gsutil cat $mod10_output_json | jq .trained_genotype_depth_depth_sepcutoff \
        | tr -d '"' >> $sub_dir/depth_cutoffs.list
        gsutil cat $mod10_output_json | jq .regeno_coverage_medians | tr -d '"' \
        >> $sub_dir/regeno_medians.list
      done < $batches_list

      # Prep input .json
      cat << EOF > cromshell/inputs/dfci-g2c.v1.$sub_name.inputs.json
{
    "RegenotypeCNVs.RD_depth_sepcutoffs": $( collapse_txt $sub_dir/depth_cutoffs.list ),
    "RegenotypeCNVs.batch_depth_vcfs": $( collapse_txt $sub_dir/depth_vcfs.list ),
    "RegenotypeCNVs.batches": $( collapse_txt $batches_list ),
    "RegenotypeCNVs.cohort": "dfci-g2c.v1",
    "RegenotypeCNVs.cohort_depth_vcf": "$staging_prefix/09/dfci-g2c.v1.all_batches.depth.vcf.gz",
    "RegenotypeCNVs.contig_list": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs.list",
    "RegenotypeCNVs.coveragefile_idxs": $( collapse_txt $sub_dir/bincov_idxs.list ),
    "RegenotypeCNVs.coveragefiles": $( collapse_txt $sub_dir/bincovs.list ),
    "RegenotypeCNVs.depth_vcfs": $( collapse_txt $sub_dir/depth_gt_vcfs.list ),
    "RegenotypeCNVs.medianfiles": $( collapse_txt $sub_dir/cov_medians.list ),
    "RegenotypeCNVs.n_RdTest_bins": 100000,
    "RegenotypeCNVs.n_per_split": 5000,
    "RegenotypeCNVs.regeno_coverage_medians": $( collapse_txt $sub_dir/regeno_medians.list ),
    "RegenotypeCNVs.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "RegenotypeCNVs.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-11-15-v1.0-488d7cb0"
}
EOF
      ;;

    #############
    # MODULE 12 #
    #############
    12)
      wdl="code/wdl/gatk-sv/CombineBatches.wdl"

      # Get URIs per batch
      while read bid; do                
        mod10_output_json=$staging_prefix/10/$bid/$bid.gatksv_module_10.outputs.json
        gsutil cat $mod10_output_json | jq .genotyped_pesr_vcf | tr -d '"' \
        >> $sub_dir/pesr_vcfs.list
        gsutil cat $mod10_output_json | jq .sr_background_fail | tr -d '"' \
        >> $sub_dir/sr_bg_fails.list
        gsutil cat $mod10_output_json | jq .sr_bothside_pass | tr -d '"' \
        >> $sub_dir/bothsides_pass.list

        gsutil ls "$staging_prefix/11/$bid.depth.regeno_final.vcf.gz" \
        >> $sub_dir/depth_vcfs.list
      done < $batches_list

      # Prep input .json template
      cat << EOF > $sub_dir/dfci-g2c.v1.$sub_name.inputs.template.json
{
    "CombineBatches.batches": $( collapse_txt $batches_list ),
    "CombineBatches.cohort_name": "dfci-g2c.v1",
    "CombineBatches.contig_list": "gs://dfci-g2c-refs/hg38/contig_fais/\$CONTIG.fai",
    "CombineBatches.depth_exclude_list": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/depth_blacklist.sorted.bed.gz",
    "CombineBatches.depth_vcfs": $( collapse_txt $sub_dir/depth_vcfs.list ),
    "CombineBatches.empty_file": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/empty.file",
    "CombineBatches.localize_shard_size": 20000,
    "CombineBatches.min_sr_background_fail_batches": 0.5,
    "CombineBatches.pe_exclude_list": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz",
    "CombineBatches.pesr_vcfs": $( collapse_txt $sub_dir/pesr_vcfs.list ),
    "CombineBatches.raw_sr_background_fail_files": $( collapse_txt $sub_dir/sr_bg_fails.list ),
    "CombineBatches.raw_sr_bothside_pass_files": $( collapse_txt $sub_dir/bothsides_pass.list ),
    "CombineBatches.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "CombineBatches.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-11-15-v1.0-488d7cb0",
    "CombineBatches.use_hail": false,
    "CombineBatches.runtime_override_sort_merged_vcf_cluster": {"disk_gb" : 50},
    "CombineBatches.runtime_override_svtk_vcf_cluster" : {"mem_gb": 15, "cpu_cores": 4, "boot_disk_gb": 15, "disk_gb": 20}
}
EOF
      ;;

    #############
    # MODULE 13 #
    #############
    13)
      wdl="code/wdl/gatk-sv/ResolveComplexVariants.wdl"

      # Collect batch-specific URIs and other info
      while read bid; do
        
        mod04_output_json=$staging_prefix/04/$bid/$bid.gatksv_module_04.outputs.json
        gsutil cat $mod04_output_json | jq .merged_PE | tr -d '"' \
        >> $sub_dir/pe_files.list

        mod07_output_json=$staging_prefix/07/$bid/$bid.gatksv_module_07.outputs.json
        gsutil cat $mod07_output_json | jq .cutoffs | tr -d '"' \
        >> $sub_dir/cutoffs.list

      done < $batches_list

      # Collect per-chromosome outputs from module 11
      for contig in $( seq 1 22 ) X Y; do
        mod12_output_json=$staging_prefix/12/chr$contig/12-CombineBatches.chr$contig.outputs.json
        gsutil cat $mod12_output_json | jq .combined_vcfs | tr -d '"' >> $sub_dir/combined_vcfs.list
        gsutil cat $mod12_output_json | jq .cluster_bothside_pass_lists \
        | tr -d '"' >> $sub_dir/bothsides_pass.list
        gsutil cat $mod12_output_json | jq .cluster_background_fail_lists \
        | tr -d '"' >> $sub_dir/background_fail.list
      done

      # Prep input .json
      cat << EOF > cromshell/inputs/dfci-g2c.v1.$sub_name.inputs.json
{
    "ResolveComplexVariants.cluster_background_fail_lists": $( collapse_txt $sub_dir/background_fail.list ),
    "ResolveComplexVariants.cluster_bothside_pass_lists": $( collapse_txt $sub_dir/bothsides_pass.list ),
    "ResolveComplexVariants.cluster_vcfs": $( collapse_txt $sub_dir/combined_vcfs.list ),
    "ResolveComplexVariants.cohort_name": "dfci-g2c.v1",
    "ResolveComplexVariants.contig_list": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs.list",
    "ResolveComplexVariants.cytobands": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/cytobands_hg38.bed.gz",
    "ResolveComplexVariants.disc_files": $( collapse_txt $sub_dir/pe_files.list ),
    "ResolveComplexVariants.max_shard_size": 500,
    "ResolveComplexVariants.mei_bed": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz",
    "ResolveComplexVariants.pe_exclude_list": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz",
    "ResolveComplexVariants.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
    "ResolveComplexVariants.rf_cutoff_files": $( collapse_txt $sub_dir/cutoffs.list ),
    "ResolveComplexVariants.runtime_override_breakpoint_overlap_filter": { "boot_disk_gb" : 20, "disk_gb" : 1000, "mem_gb" : 31.5, "preemptible_tries" : 1 },
    "ResolveComplexVariants.runtime_override_concat_resolved_per_shard": { "boot_disk_gb" : 20, "disk_gb" : 150, "mem_gb" : 15.5, "preemptible_tries" : 1, "cpu_cores" : 4 },
    "ResolveComplexVariants.runtime_override_integrate_resolved_vcfs": { "boot_disk_gb" : 20, "disk_gb" : 500 },
    "ResolveComplexVariants.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "ResolveComplexVariants.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-11-15-v1.0-488d7cb0"
}
EOF
      ;;

    #############
    # MODULE 14 #
    #############
    14)
      wdl="code/wdl/gatk-sv/GenotypeComplexVariants.wdl"

      # Get URIs per batch
      while read bid; do
        mod04_output_json=$staging_prefix/04/$bid/$bid.gatksv_module_04.outputs.json
        gsutil cat $mod04_output_json | jq .merged_bincov | tr -d '"' \
        >> $sub_dir/bincovs.list

        mod10_output_json=$staging_prefix/10/$bid/$bid.gatksv_module_10.outputs.json
        gsutil cat $mod10_output_json | jq .trained_genotype_depth_depth_sepcutoff \
        | tr -d '"' >> $sub_dir/depth_cutoffs.list

        mod14A_output_json=$staging_prefix/14A/$bid/$bid.gatksv_module_14A.outputs.json
        gsutil cat $mod14A_output_json | jq .filtered_median_cov | tr -d '"' \
        >> $sub_dir/median_covs.list

        echo "$staging_prefix/11/$bid.depth.regeno_final.vcf.gz" \
        >> $sub_dir/depth_regeno_vcfs.list
      done < $batches_list

      # Get arrays of inputs per chromosome
      for k in $( seq 1 22 ) X Y; do
        echo -e "$staging_prefix/13/dfci-g2c.v1.reshard_vcf.chr$k.resharded.vcf.gz.tbi" \
        >> $sub_dir/cpx_vcf_idxs.list
        echo -e "$staging_prefix/13/dfci-g2c.v1.reshard_vcf.chr$k.resharded.vcf.gz" \
        >> $sub_dir/cpx_vcfs.list
      done

      # Prep input .json
      cat << EOF > cromshell/inputs/dfci-g2c.v1.$sub_name.inputs.json
{
    "GenotypeComplexVariants.batches": $( collapse_txt $batches_list ),
    "GenotypeComplexVariants.bin_exclude": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/bin_exclude.hg38.gatkcov.bed.gz",
    "GenotypeComplexVariants.bincov_files": $( collapse_txt $sub_dir/bincovs.list ),
    "GenotypeComplexVariants.cohort_name": "dfci-g2c.v1",
    "GenotypeComplexVariants.complex_resolve_vcf_indexes": $( collapse_txt $sub_dir/cpx_vcf_idxs.list ),
    "GenotypeComplexVariants.complex_resolve_vcfs": $( collapse_txt $sub_dir/cpx_vcfs.list ),
    "GenotypeComplexVariants.contig_list": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs.list",
    "GenotypeComplexVariants.depth_gt_rd_sep_files": $( collapse_txt $sub_dir/depth_cutoffs.list ),
    "GenotypeComplexVariants.depth_vcfs": $( collapse_txt $sub_dir/depth_regeno_vcfs.list ),
    "GenotypeComplexVariants.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
    "GenotypeComplexVariants.median_coverage_files": $( collapse_txt $sub_dir/median_covs.list ),
    "GenotypeComplexVariants.ped_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.gatksv_module14.ped",
    "GenotypeComplexVariants.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
    "GenotypeComplexVariants.runtime_override_parse_genotypes": {"disk_gb" : 275, "boot_disk_gb" : 30, "mem_gb" : 7.5, "cpu_cores" : 2},
    "GenotypeComplexVariants.runtime_override_rd_genotype": {"disk_gb" : 20, "mem_gb" : 15.5, "n_cpu" : 4, "preemptible_tries" : 1},
    "GenotypeComplexVariants.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "GenotypeComplexVariants.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-01-14-v1.0.1-88dbd052"
}
EOF
      ;;

    #############
    # MODULE 15 #
    #############
    15)
      wdl="code/wdl/gatk-sv/CleanVcf.wdl"

      # Get arrays of inputs per chromosome
      for k in $( seq 1 22 ) X Y; do
        echo -e "$staging_prefix/14/dfci-g2c.v1.chr$k.regenotyped.vcf.gz" \
        >> $sub_dir/vcfs.list
      done

      # Prep input .json
      cat << EOF > cromshell/inputs/dfci-g2c.v1.$sub_name.inputs.json
{
    "CleanVcf.HERVK_reference": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/HERVK.sorted.bed.gz",
    "CleanVcf.LINE1_reference": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/LINE1.sorted.bed.gz",
    "CleanVcf.allosome_fai": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/allosome.fai",
    "CleanVcf.chr_x": "chrX",
    "CleanVcf.chr_y": "chrY",
    "CleanVcf.clean_vcf1b_records_per_shard": 10000,
    "CleanVcf.clean_vcf5_records_per_shard": 5000,
    "CleanVcf.cohort_name": "dfci-g2c.v1",
    "CleanVcf.complex_genotype_vcfs": $( collapse_txt $sub_dir/vcfs.list ),
    "CleanVcf.complex_resolve_background_fail_list": "$staging_prefix/13/dfci-g2c.v1.all.sr_background_fail.updated3.txt",
    "CleanVcf.complex_resolve_bothside_pass_list": "$staging_prefix/13/dfci-g2c.v1.all.sr_bothside_pass.updated3.txt",
    "CleanVcf.contig_list": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/contig.fai",
    "CleanVcf.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
    "CleanVcf.max_shards_per_chrom_step1": 200,
    "CleanVcf.min_records_per_shard_step1": 5000,
    "CleanVcf.ped_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.gatksv_module14.ped",
    "CleanVcf.primary_contigs_list": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs.list",
    "CleanVcf.runtime_attr_override_build_dict_1b": {"max_retries" : 1, "preemptible_tries" : 1, "mem_gb" : 15.5, "cpu_cores" : 4, "disk_gb" : 100, "boot_disk_gb" : 30},
    "CleanVcf.runtime_override_clean_vcf_5_polish": {"max_retries" : 1, "preemptible_tries" : 1, "mem_gb" : 128, "cpu_cores" : 16, "disk_gb" : 550, "boot_disk_gb" : 30},
    "CleanVcf.runtime_override_split_vcf_to_clean": {"max_retries" : 1, "preemptible_tries" : 1, "disk_gb" : 1000, "boot_disk_gb" : 20},
    "CleanVcf.samples_per_step2_shard": 100,
    "CleanVcf.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "CleanVcf.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-01-14-v1.0.1-88dbd052"
}
EOF
      ;;

    #############
    # MODULE 16 #
    #############
    16)
      wdl="code/wdl/gatk-sv/RefineComplexVariants.wdl"

      # Get batch-specific inputs
      while read bid; do
        # Batch membership
        echo $bid >> $sub_dir/batch_names.list
        echo \
          "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/module14A_batch_sample_lists/$bid.gatksv_module14A.samples.list" \
        >> $sub_dir/batch_sample_lists.list

        # Module 04 outputs
        mod04_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/04/$bid/$bid.gatksv_module_04.outputs.json
        gsutil cat $mod04_json | jq .merged_PE | tr -d '"' \
        >> $sub_dir/pe_evidence.list
        gsutil cat $mod04_json | jq .merged_PE_index | tr -d '"' \
        >> $sub_dir/pe_evidence_idx.list
        gsutil cat $mod04_json | jq .merged_dels | tr -d '"' \
        >> $sub_dir/depth_dels.list
        gsutil cat $mod04_json | jq .merged_dups | tr -d '"' \
        >> $sub_dir/depth_dups.list
      done < $batches_list

      # Prep input .json
      cat << EOF > $sub_dir/dfci-g2c.v1.$sub_name.inputs.template.json
{
    "RefineComplexVariants.Depth_DEL_beds": $( collapse_txt $sub_dir/depth_dels.list ),
    "RefineComplexVariants.Depth_DUP_beds": $( collapse_txt $sub_dir/depth_dups.list ),
    "RefineComplexVariants.GenerateCpxReviewScript.script": "$MAIN_WORKSPACE_BUCKET/code/scripts/reformat_CPX_bed_and_generate_script.g2c.py",
    "RefineComplexVariants.PE_metrics": $( collapse_txt $sub_dir/pe_evidence.list ),
    "RefineComplexVariants.PE_metrics_indexes": $( collapse_txt $sub_dir/pe_evidence_idx.list ),
    "RefineComplexVariants.batch_name_list": $( collapse_txt $sub_dir/batch_names.list ),
    "RefineComplexVariants.batch_sample_lists": $( collapse_txt $sub_dir/batch_sample_lists.list ),
    "RefineComplexVariants.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
    "RefineComplexVariants.n_per_split": 15000,
    "RefineComplexVariants.prefix": "dfci-g2c.v1.\$CONTIG",
    "RefineComplexVariants.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "RefineComplexVariants.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-01-14-v1.0.1-88dbd052",
    "RefineComplexVariants.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/15/dfci-g2c.v1.\$CONTIG.final_format.vcf.gz"
}
EOF
      ;;

    #############
    # MODULE 17 #
    #############
    17)
      wdl="code/wdl/gatk-sv/JoinRawCalls.wdl"

      # For speed, since so many jq calls are required, we can copy all
      # output json files from module 05 into the submission directory
      gsutil -m cp \
        $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/05C/*/*.gatksv_module_05C.outputs.json \
        $sub_dir/

      # Get URIs per batch
      while read bid; do

        mod05C_output_json=$sub_dir/$bid.gatksv_module_05C.outputs.json

        for alg in depth manta wham melt; do
          jq .clustered_${alg}_vcf $mod05C_output_json | tr -d '"' \
          >> $sub_dir/${alg}_vcfs.list
          jq .clustered_${alg}_vcf_index $mod05C_output_json | tr -d '"' \
          >> $sub_dir/${alg}_vcf_idxs.list
        done

      done < $batches_list

      # Prep input .json template
      cat << EOF > cromshell/inputs/dfci-g2c.v1.$sub_name.inputs.json
{
    "JoinRawCalls.FormatVcfForGatk.formatter_args": "--fix-end",
    "JoinRawCalls.clustered_depth_vcfs": $( collapse_txt $sub_dir/depth_vcfs.list ),
    "JoinRawCalls.clustered_depth_vcf_indexes": $( collapse_txt $sub_dir/depth_vcf_idxs.list ),
    "JoinRawCalls.clustered_manta_vcfs":  $( collapse_txt $sub_dir/manta_vcfs.list ),
    "JoinRawCalls.clustered_manta_vcf_indexes": $( collapse_txt $sub_dir/manta_vcf_idxs.list ),
    "JoinRawCalls.clustered_melt_vcfs":  $( collapse_txt $sub_dir/melt_vcfs.list ),
    "JoinRawCalls.clustered_melt_vcf_indexes": $( collapse_txt $sub_dir/melt_vcf_idxs.list ),
    "JoinRawCalls.clustered_wham_vcfs":  $( collapse_txt $sub_dir/wham_vcfs.list ),
    "JoinRawCalls.clustered_wham_vcf_indexes": $( collapse_txt $sub_dir/wham_vcf_idxs.list ),
    "JoinRawCalls.contig_list": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs.list",
    "JoinRawCalls.gatk_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2024-12-05-4.6.1.0-6-gfc248dfc1-NIGHTLY-SNAPSHOT",
    "JoinRawCalls.ped_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
    "JoinRawCalls.prefix": "dfci-g2c.v1",
    "JoinRawCalls.reference_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
    "JoinRawCalls.reference_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
    "JoinRawCalls.reference_fasta_fai": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
    "JoinRawCalls.runtime_attr_svcluster": { "disk_gb" : 250, "mem_gb" : 128, cpu_cores : 32, "boot_disk_gb" : 20 },
    "JoinRawCalls.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "JoinRawCalls.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-01-14-v1.0.1-88dbd052"
}
EOF
      ;;

    #############
    # MODULE 18 #
    #############
    18)
      wdl="code/wdl/gatk-sv/SVConcordance.wdl"

      # Prep input .json
      cat << EOF > $sub_dir/dfci-g2c.v1.$sub_name.inputs.template.json
{
    "SVConcordance.contig_list": "gs://dfci-g2c-refs/hg38/contig_lists/\$CONTIG.list",
    "SVConcordance.eval_vcf": "$staging_prefix/PosthocHardFilterPart2/\$CONTIG/HardFilterPart2/dfci-g2c.v1.\$CONTIG.cpx_refined.posthoc_filtered.posthoc_filtered.vcf.gz",
    "SVConcordance.gatk_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2024-12-05-4.6.1.0-6-gfc248dfc1-NIGHTLY-SNAPSHOT",
    "SVConcordance.output_prefix": "dfci-g2c.v1.\$CONTIG",
    "SVConcordance.reference_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
    "SVConcordance.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "SVConcordance.truth_vcf": "$staging_prefix/17/dfci-g2c.v1.join_raw_calls.vcf.gz"
}
EOF
      ;;

    #############
    # MODULE 19 #
    #############
    19)
      wdl="code/wdl/gatk-sv/FilterGenotypes.wdl"

      # Prep input .json
      cat << EOF > $sub_dir/dfci-g2c.v1.$sub_name.inputs.template.json
{
    "FilterGenotypes.gatk_docker": "us.gcr.io/broad-dsde-methods/markw/gatk:mw-tb-form-sv-filter-training-data-899360a",
    "FilterGenotypes.genome_tracks": ["gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-RepeatMasker.bed.gz",
                                      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-Segmental-Dups.bed.gz",
                                      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-Simple-Repeats.bed.gz",
                                      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38_umap_s100.bed.gz",
                                      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38_umap_s24.bed.gz"],
    "FilterGenotypes.gq_recalibrator_model_file": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/gatk-sv-recalibrator.aou_phase_1.v1.model",
    "FilterGenotypes.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
    "FilterGenotypes.no_call_rate_cutoff": 1,
    "FilterGenotypes.output_prefix": "dfci-g2c.v1.\$CONTIG",
    "FilterGenotypes.ped_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
    "FilterGenotypes.ploidy_table": "$staging_prefix/17/dfci-g2c.v1.ploidy.tsv",
    "FilterGenotypes.primary_contigs_fai": "gs://dfci-g2c-refs/hg38/contig_fais/\$CONTIG.fai",
    "FilterGenotypes.recalibrate_gq_args": ["--keep-homvar false","--keep-homref true","--keep-multiallelic true","--skip-genotype-filtering true","--min-samples-to-estimate-allele-frequency -1"],
    "FilterGenotypes.run_qc": false,
    "FilterGenotypes.runtime_override_plot_qc_per_family": {"mem_gb" : 15, "disk_gb" : 100},
    "FilterGenotypes.sl_filter_args": "",
    "FilterGenotypes.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "FilterGenotypes.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-01-14-v1.0.1-88dbd052",
    "FilterGenotypes.vcf": "$staging_prefix/18/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.concordance.vcf.gz"
}
EOF
      ;;

    *)
      echo "Module number $module_idx not recognized by submit_cohort_module. Exiting."
      return 2
      ;;

  esac

  # Submit job and add job ID to list of jobs for this module
  case $module_idx in
    # Submission for modules 12, 16, and 18-19 are handled differently due to them being contig-sharded
    12|16|18|19)
      code/scripts/manage_chromshards.py \
        --wdl $wdl \
        --input-json-template $sub_dir/dfci-g2c.v1.$sub_name.inputs.template.json \
        --dependencies-zip gatksv.dependencies.zip \
        --staging-bucket $staging_prefix/$module_idx \
        --name $sub_name \
        --status-tsv cromshell/progress/dfci-g2c.v1.$sub_name.progress.tsv \
        --workflow-id-log-prefix "dfci-g2c.v1" \
        --outer-gate 45 \
        --max-attempts $max_attempts
      ;;

    *)
      cmd="cromshell --no_turtle -t 120 -mc submit"
      cmd="$cmd --options-json code/refs/json/aou.cromwell_options.default.json"
      cmd="$cmd --dependencies-zip gatksv.dependencies.zip"
      cmd="$cmd $wdl cromshell/inputs/dfci-g2c.v1.$sub_name.inputs.json"
      eval $cmd | jq .id | tr -d '"' \
      >> cromshell/job_ids/dfci-g2c.v1.$sub_name.job_ids.list
      ;;
  esac
}


#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Main variant filtering after callset integration

# Note that this code is designed to run on the All of Us Workbench v2.0 (Verily Pre)
# It may not have reverse compatability with the legacy All of Us Workbench v1.0


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://rw-migration-aou-rw-84a0039b
gcloud storage cp $MAIN_WORKSPACE_BUCKET/code/scripts/configure_verily_vm.sh ./ && \
. configure_verily_vm.sh && \
rm configure_verily_vm.sh

# Set up local directory structure
for dir in cromwell cromwell/inputs cromwell/submissions; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done


##################################################
# Curate sample-level features for all filtering #
##################################################

# TODO: implement this


##############################################
# Curate SR/SD mask for filtering annotation #
##############################################

# Reaffirm staging directory
staging_dir=staging/sd_sr_curation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Download hg38 RepMask from UCSC and filter to relevant subset
# TODO: FINISH THIS
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz



###################################################################
# Curate variant- and genotype-level features for SV GT filtering #
###################################################################

# BELOW IS DEV: submit as a single job to familiarize with wb workflow syntax

# Must run to register workflow in workspace (need to re-register every time code changes)
wdl_version=$( gsutil cat $MAIN_WORKSPACE_BUCKET/code/refs/wdl.version_info.txt \
               | awk '{ if ($1=="pancan_germline_wgs") print $3 }' )
workflow_base="CollectGTFilterFeatures"
workflow_name="${workflow_base}_${wdl_version}"
if [ $( wb workflow list | fgrep $workflow_name | wc -l ) -lt 1 ]; then
  wb workflow create \
    --bucket-id main_g2c_bucket \
    --path code/wdl/pancan_germline_wgs/CollectGTFilterFeatures.wdl \
    --workflow "$workflow_name" \
    --display-name "$workflow_name" \
    --workflow-type WDL
fi

# Create inputs and upload to workspace bucket
input_json="$workflow_name.inputs.$( date +"%Y%m%d_%H%M%S" ).json"
cat << EOF > "cromwell/inputs/$input_json"
{
  "CollectGTFilterFeatures.g2c_analysis_docker": "vanallenlab/g2c_analysis:84838e6",
  "CollectGTFilterFeatures.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/indel_sv_integration/chr22/ExtractLargeSvs/shard-2/dfci-g2c.v1.sv.chr22.3.sorted.vcf.gz",
  "CollectGTFilterFeatures.vcf_idx": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/indel_sv_integration/chr22/ExtractLargeSvs/shard-2/dfci-g2c.v1.sv.chr22.3.sorted.vcf.gz.tbi"
}
EOF
gsutil cp "cromwell/inputs/$input_json" $WORKSPACE_BUCKET/cromwell-inputs/$workflow_base/

# Test submission for a single chr22 SV VCF
wb workflow job run \
  --workflow "$workflow_name" \
  --inputs-uri "$WORKSPACE_BUCKET/cromwell-inputs/$workflow_base/$input_json" \
  --output-bucket-id "$( echo $WORKSPACE_BUCKET | sed 's/^gs\:\/\///g'  )" \
  --output-path cromwell-execution \
  --read-from-cache \
  --write-to-cache \
  --format JSON \
| python -m json.tool \
> "cromwell/submissions/$workflow_name.submission.json"

# Check progress
wid=$( jq .runId "cromwell/submissions/$workflow_name.submission.json" | tr -d '"' )
wb workflow job describe --job-id $wid --format JSON \
| python -m json.tool

# TODO: implement this



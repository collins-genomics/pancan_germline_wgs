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


################################################
# Curate sample-level covariates for filtering #
################################################

# BELOW IS DEV: submit as a single job to familiarize with wb workflow syntax
# Must run to register workflow in workspace (need to rerun every time code changes)
wb workflow create \
  --bucket-id main_g2c_bucket \
  --path code/wdl/pancan_germline_wgs/CollectGTFilterFeatures.wdl \
  --workflow CollectGTFilterFeatures_e8b4d27 \
  --display-name CollectGTFilterFeatures_e8b4d27 \
  --workflow-type WDL

# Create inputs and upload to workspace bucket
cat << EOF > ~/scratch/CollectGTFilterFeatures.test.inputs.json
{
  "CollectGTFilterFeatures.g2c_analysis_docker": "vanallenlab/g2c_analysis:84838e6",
  "CollectGTFilterFeatures.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/indel_sv_integration/chr22/ExtractLargeSvs/shard-2/dfci-g2c.v1.sv.chr22.3.sorted.vcf.gz",
  "CollectGTFilterFeatures.vcf_idx": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/indel_sv_integration/chr22/ExtractLargeSvs/shard-2/dfci-g2c.v1.sv.chr22.3.sorted.vcf.gz.tbi"
}
EOF

# Test submission for a single chr22 SV VCF
wb workflow job run \
  --workflow CollectGTFilterFeatures_e8b4d27 \
  --inputs-uri TBD \
  --output-bucket-id $( echo $WORKSPACE_BUCKET | sed 's/^gs\:\/\///g'  ) \
  --output-path cromwell-output/ \
  --read-from-cache \
  --write-to-cache


# TODO: implement this


##############################################################
# Compute sample-level global no-call rate per variant class #
##############################################################

# TODO: implement this


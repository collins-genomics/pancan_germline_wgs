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
for dir in cromwell cromwell/inputs cromwell/submissions staging; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done


##################################################
# Curate sample-level features for all filtering #
##################################################

# TODO: implement this


###############################################
# Curate UCSC tracks for filtering annotation #
###############################################

# NOTE: this chunk must be run outside of AoU Verily workbench (can be local, files are small)

# Reaffirm staging directory
staging_dir=~/staging/ucsc_curation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Download hg38 RepMask from UCSC and filter to relevant subset
wget -O - \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz \
| gunzip -c \
| awk -v FS="\t" -v OFS="\t" \
  '{ if ($12 ~ /Simple_repeat|Satellite|Low_complexity/) print $6, $7, $8 }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bgzip -c \
> $staging_dir/hg38.simple_repeats.bed.gz
tabix -p bed -f $staging_dir/hg38.simple_repeats.bed.gz

# Download hg38 segdups from UCSC and curate
wget -O - \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz \
| gunzip -c \
| awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4 }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bgzip -c \
> $staging_dir/hg38.segdups.bed.gz
tabix -p bed -f $staging_dir/hg38.segdups.bed.gz

# Combine SR and SD tracks
zcat \
  $staging_dir/hg38.simple_repeats.bed.gz \
  $staging_dir/hg38.segdups.bed.gz \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bgzip -c \
> $staging_dir/hg38.sr_sd.bed.gz
tabix -p bed -f $staging_dir/hg38.sr_sd.bed.gz

# Also download and stage 100mer UMap bigwig for convenience
wget \
  -O $staging_dir/hg38.100mer_multiUmap.bw \
  http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw

# Copy all repeat tracks to ref bucket for future reference
gsutil -m cp \
  $staging_dir/hg38.*.bed.gz* \
  $staging_dir/hg38.100mer_multiUmap.bw \
  gs://dfci-g2c-refs/ucsc/hg38/

# Clean up garbage
rm -rf $staging_dir


###################################################################
# Curate variant- and genotype-level features for SV GT filtering #
###################################################################

# BELOW IS DEV: manual submit for chr22 to familiarize with wb workflow syntax

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

# Reaffirm staging directory
staging_dir=~/staging/sv_gt_feature_curation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Dev: set contig manually (this can be looped once automated)
contig=chr22

# Create VCF input array
# NOTE: because this was the first task after migrating from AoU RW v1.0 to Verily,
# we need to update the URI prefixes for all VCFs to the new Verily buckets
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/indel_sv_integration/$contig/UnifyGatkCallsets.$contig.outputs.json \
| jq .\"UnifyGatkCallsets.cleaned_sv_vcfs\" \
| fgrep "gs://" | tr -d '",' \
| awk -v FS="42de3578da45" -v prefix="$MAIN_WORKSPACE_BUCKET" -v OFS="\t" \
  '{ print prefix$2, prefix$2".tbi" }' \
> $staging_dir/$workflow_name.vcf_inputs.$contig.tsv
gsutil -m cp \
  $staging_dir/$workflow_name.vcf_inputs.$contig.tsv \
  $MAIN_WORKSPACE_BUCKET/data/sv_gt_filtering/

# Create inputs and upload to workspace bucket
input_json="$workflow_name.$contig.inputs.$( date +"%Y%m%d_%H%M%S" ).json"
cat << EOF > "cromwell/inputs/$input_json"
{
  "CollectGTFilterFeatures.bed_features": ["gs://dfci-g2c-refs/giab/$contig/giab.hg38.broad_callable.hard.$contig.bed.gz",
                                           "gs://dfci-g2c-refs/ucsc/hg38/hg38.segdups.bed.gz",
                                           "gs://dfci-g2c-refs/ucsc/hg38/hg38.simple_repeats.bed.gz"],
  "CollectGTFilterFeatures.bigwig_features": ["gs://dfci-g2c-refs/ucsc/hg38/hg38.100mer_multiUmap.bw"],
  "CollectGTFilterFeatures.bigwig_feature_names": ["umap"],
  "CollectGTFilterFeatures.bed_feature_names": ["giab_hard", "segdup", "simrep"],
  "CollectGTFilterFeatures.g2c_analysis_docker": "vanallenlab/g2c_analysis:2d676dc",
  "CollectGTFilterFeatures.linux_docker": "ubuntu:plucky-20251001",
  "CollectGTFilterFeatures.output_prefix": "dfci-g2c.v1.$contig.sv",
  "CollectGTFilterFeatures.vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/data/sv_gt_filtering/$workflow_name.vcf_inputs.$contig.tsv"
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
> "cromwell/submissions/$workflow_name.$contig.submission.json"

# Check progress
wid=$( jq .runId "cromwell/submissions/$workflow_name.$contig.submission.json" | tr -d '"' )
wb workflow job describe --job-id $wid --format JSON \
| python -m json.tool | jq .status

# TODO: finish implementing this



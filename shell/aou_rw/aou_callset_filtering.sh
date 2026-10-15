#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Main variant filtering after callset integration

# Note that this code is designed to run on the All of Us Workbench v2.0 (Verily Pre)
# It may not have reverse compatability with the legacy All of Us Workbench v1.0

# The below code is responsible for generating two parallel filtered versions of the G2C callset:
# 1. "Alpha" (v0.9): direct application of existing filter methods
# 2. "Main" (v1.0): development & application of a brand new universal filtering method
# These two versions ("alpha" vs "main") are referenced shorthand below


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
for dir in staging; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy legacy mingq code
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/misc/legacy_mingq_wdl code/wdl/gatk-sv/

# Create dependencies .zip for GATK-SV module submissions
cd code/wdl/gatk-sv && \
zip -r gatksv.dependencies.zip . && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for minGQ module submissions
cd code/wdl/gatk-sv/legacy_mingq_wdl && \
zip -r mingq.dependencies.zip . && \
mv mingq.dependencies.zip ~/ && \
cd ~

# Launch Cromwell in server mode on the same Verily VM
# Note that this must be executed from a separate terminal after 
# configuring the environment per the above:
# code/scripts/launch_cromwell_server.sh


# ####################################################################
# # [Alpha] Apply pre-trained high-sensitivity FilterGenotypes model #
# ####################################################################

# # Make staging directory
# staging_dir=staging/FilterGenotypes
# if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for GATK-SV module 19
cat << EOF > $staging_dir/FilterGenotypes.inputs.template.json
{
  "FilterGenotypes.FilterVcf.script": "$MAIN_WORKSPACE_BUCKET/code/scripts/apply_sl_filter.g2c_imputed.py",
  "FilterGenotypes.gatk_docker": "us.gcr.io/broad-dsde-methods/markw/gatk:mw-tb-form-sv-filter-training-data-899360a",
  "FilterGenotypes.genome_tracks": ["gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-RepeatMasker.bed.gz",
                                    "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-Segmental-Dups.bed.gz",
                                    "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-Simple-Repeats.bed.gz",
                                    "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38_umap_s100.bed.gz",
                                    "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38_umap_s24.bed.gz"],
  "FilterGenotypes.gq_recalibrator_model_file": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/gatk-sv-recalibrator.aou_phase_1.v1.model",
  "FilterGenotypes.linux_docker": "ubuntu:plucky-20251001",
  "FilterGenotypes.no_call_rate_cutoff": 1,
  "FilterGenotypes.output_prefix": ""dfci-g2c.v0.9.\$CONTIG",
  "FilterGenotypes.ped_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "FilterGenotypes.ploidy_table": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/17/dfci-g2c.v1.ploidy.tsv",
  "FilterGenotypes.primary_contigs_fai": "gs://dfci-g2c-refs/hg38/contig_fais/\$CONTIG.fai",
  "FilterGenotypes.recalibrate_gq_args": ["--keep-homvar false","--keep-homref true","--keep-multiallelic true","--skip-genotype-filtering true","--min-samples-to-estimate-allele-frequency -1"],
  "FilterGenotypes.run_qc": false,
  "FilterGenotypes.runtime_override_plot_qc_per_family": {"mem_gb" : 15, "disk_gb" : 100},
  "FilterGenotypes.sl_filter_args": "--small-del-threshold -53 --medium-del-threshold -12 --small-dup-threshold -105 --medium-dup-threshold -81 --ins-threshold -97",
  "FilterGenotypes.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "FilterGenotypes.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-06-27-v1.0.4-63e6c81e",
  "FilterGenotypes.vcf": "$WORKSPACE_BUCKET/qc/header_typo_fix/\$CONTIG/FixTypo/dfci-ufc.gatksv.v2.typo_fixed.vcf.gz"
}
EOF

# # Run module 19 on all UFC SV VCFs
# code/scripts/manage_chromshards.py \
#   --wdl code/wdl/gatk-sv/FilterGenotypes.wdl \
#   --input-json-template $staging_dir/FilterGenotypes.inputs.template.json \
#   --staging-bucket $WORKSPACE_BUCKET/data/module19 \
#   --dependencies-zip gatksv.dependencies.zip \
#   --status-tsv cromshell/progress/FilterGenotypes.progress.tsv \
#   --workflow-id-log-prefix "dfci-ufc.v2" \
#   --outer-gate 30 \
#   --submission-gate 0.1 \
#   --max-attempts 2


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
> $staging_dir/CollectGTFilterFeatures.vcf_inputs.$contig.tsv
gsutil -m cp \
  $staging_dir/CollectGTFilterFeatures.vcf_inputs.$contig.tsv \
  $MAIN_WORKSPACE_BUCKET/data/sv_gt_filtering/

# # Create inputs and upload to workspace bucket
# input_json="CollectGTFilterFeatures.$contig.inputs.json"
# cat << EOF > "cromwell/inputs/$input_json"
# {
#   "CollectGTFilterFeatures.bed_features": ["gs://dfci-g2c-refs/giab/$contig/giab.hg38.broad_callable.hard.$contig.bed.gz",
#                                            "gs://dfci-g2c-refs/ucsc/hg38/hg38.segdups.bed.gz",
#                                            "gs://dfci-g2c-refs/ucsc/hg38/hg38.simple_repeats.bed.gz"],
#   "CollectGTFilterFeatures.bigwig_features": ["gs://dfci-g2c-refs/ucsc/hg38/hg38.100mer_multiUmap.bw"],
#   "CollectGTFilterFeatures.bigwig_feature_names": ["umap"],
#   "CollectGTFilterFeatures.bed_feature_names": ["giab_hard", "segdup", "simrep"],
#   "CollectGTFilterFeatures.g2c_analysis_docker": "vanallenlab/g2c_analysis:77bb125",
#   "CollectGTFilterFeatures.linux_docker": "ubuntu:plucky-20251001",
#   "CollectGTFilterFeatures.output_prefix": "dfci-g2c.v1.$contig.sv",
#   "CollectGTFilterFeatures.vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/data/sv_gt_filtering/CollectGTFilterFeatures.vcf_inputs.$contig.tsv"
# }
# EOF
# gsutil cp "cromwell/inputs/$input_json" $WORKSPACE_BUCKET/cromwell-inputs/$workflow_base/

# # Test submission for a single chr22 SV VCF
# wb workflow job run \
#   --workflow "CollectGTFilterFeatures" \
#   --inputs-uri "$WORKSPACE_BUCKET/cromwell-inputs/$workflow_base/$input_json" \
#   --output-bucket-id "$( echo $WORKSPACE_BUCKET | sed 's/^gs\:\/\///g'  )" \
#   --output-path cromwell-execution \
#   --read-from-cache \
#   --write-to-cache \
#   --format JSON \
# | python -m json.tool \
# > "cromwell/submissions/CollectGTFilterFeatures.$contig.submission.json"

# Write template input .json for SV GT refinement
cat << EOF > $staging_dir/CollectSVGTFilterFeatures.inputs.template.json
{
  "CollectGTFilterFeatures.bed_features": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz",
                                           "gs://dfci-g2c-refs/ucsc/hg38/hg38.segdups.bed.gz",
                                           "gs://dfci-g2c-refs/ucsc/hg38/hg38.simple_repeats.bed.gz"],
  "CollectGTFilterFeatures.bigwig_features": ["gs://dfci-g2c-refs/ucsc/hg38/hg38.100mer_multiUmap.bw"],
  "CollectGTFilterFeatures.bigwig_feature_names": ["umap"],
  "CollectGTFilterFeatures.bed_feature_names": ["giab_hard", "segdup", "simrep"],
  "CollectGTFilterFeatures.g2c_analysis_docker": "vanallenlab/g2c_analysis:77bb125",
  "CollectGTFilterFeatures.linux_docker": "ubuntu:plucky-20251001",
  "CollectGTFilterFeatures.output_prefix": "dfci-g2c.v1.\$CONTIG.sv",
  "CollectGTFilterFeatures.vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/data/sv_gt_filtering/CollectGTFilterFeatures.vcf_inputs.\$CONTIG.tsv"
}
EOF

# Submit, monitor, and stage/cleanup SV GT refinement
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/CollectGTFilterFeatures.wdl \
  --input-json-template $staging_dir/CollectSVGTFilterFeatures.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/CollectSVGTFilterFeatures \
  --name CollectSVGTFilterFeatures \
  --contig-list <( fgrep chr22 contig_lists/dfci-g2c.v1.contigs.$WN.list ) \
  --status-tsv cromshell/progress/dfci-g2c.v1.CollectSVGTFilterFeatures.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1"



# # Check progress
# wid=$( jq .runId "cromwell/submissions/CollectGTFilterFeatures.$contig.submission.json" | tr -d '"' )
# monitor_workflow $wid
# wb workflow job describe --job-id=$wid --format=JSON | python -m json.tool

# TODO: finish implementing this



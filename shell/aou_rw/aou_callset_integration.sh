#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Variant class integration after initial callset QC

# Note that this code is designed to be run inside the AoU Researcher Workbench


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in cromshell cromshell/inputs cromshell/inputs/templates \
           cromshell/job_ids cromshell/progress staging; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}
find code/ -name "*.R" | xargs -I {} chmod a+x {}
find code/ -name "*.sh" | xargs -I {} chmod a+x {}

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Create dependencies .zip for generic G2C workflow submissions
cd code/wdl/pancan_germline_wgs && \
zip -r g2c.dependencies.zip . && \
mv g2c.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for QC workflow submissions
cd code/wdl/pancan_germline_wgs/vcf-qc && \
zip qc.dependencies.zip *.wdl && \
mv qc.dependencies.zip ~/ && \
cd ~

# Ensure Cromwell/Cromshell are configured
code/scripts/setup_cromshell.py

# Install necessary packages
. code/refs/install_packages.sh python R

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Download workspace-specific contig lists
gsutil cp -r \
  gs://dfci-g2c-refs/hg38/contig_lists \
  ./


#################################################
# Refine common SV genotypes with flanking SNVs #
#################################################

# Reaffirm staging directory
staging_dir=staging/sv_gt_cleanup
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Build SNV exclusion list
gsutil -m cat \
  gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/chr*/gnomad.v4.1.chr*.sv.sites.bed.gz \
| gunzip -c \
| awk -v FS="\t" -v OFS="\t" '{ if ($6 ~ /DEL|DUP|CNV/ && $9>=0.05) print $1, $2, $3 }' \
| cat - \
  <( gsutil -m cat \
       gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/hg38.SD_gaps_Cen_Tel_Heter_Satellite_lumpy.blacklist.sorted.merged.bed.gz \
     | gunzip -c ) \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bgzip -c \
> $staging_dir/dfci-g2c.v1.sv_regenotyping.snv_mask.bed.gz
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.sv_regenotyping.snv_mask.bed.gz \
  $MAIN_WORKSPACE_BUCKET/data/sv_regenotyping/

# Prepare sample covariates for refining imputation
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/dfci-g2c.intake_qc.all.tsv.gz \
  $staging_dir/
code/scripts/prep_covariates_for_sv_regenotyping.R \
  --qc-tsv $staging_dir/dfci-g2c.intake_qc.all.tsv.gz \
  --out-tsv $staging_dir/dfci-g2c.v1.sv_imputation_covariates.tsv
gzip -f $staging_dir/dfci-g2c.v1.sv_imputation_covariates.tsv
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.sv_imputation_covariates.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/data/sv_regenotyping/

# Make lists of all SNV VCFs and indexes for each contig
while read contig; do
  gsutil -m cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart2/$contig/PosthocCleanupPart2.$contig.outputs.json \
  | jq ".filtered_vcfs" \
  | tr -d "[]" | sed 's/,/\n/g' | fgrep "vcf.gz" | tr -d '"' | sort -V \
  | awk -v OFS="\t" '{ print $1, $1".tbi" }' \
  > $staging_dir/dfci-g2c.v1.sv_regenotyping.snv_vcf_info.$contig.tsv
done < contig_lists/dfci-g2c.v1.contigs.$WN.list
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.sv_regenotyping.snv_vcf_info.*.tsv \
  $MAIN_WORKSPACE_BUCKET/data/sv_regenotyping/

# Write template input .json for SV GT refinement
cat << EOF > $staging_dir/RefineSvGenotypesWithSnvs.inputs.template.json
{
  "RefineSvGenotypesWithSnvs.QuerySnvs.n_preemptible": 0,
  "RefineSvGenotypesWithSnvs.tmp_dev_docker": "vanallenlab/g2c_analysis:dd6cccc",
  "RefineSvGenotypesWithSnvs.g2c_analysis_docker": "vanallenlab/g2c_analysis:dd6cccc",
  "RefineSvGenotypesWithSnvs.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "RefineSvGenotypesWithSnvs.linux_docker": "ubuntu:plucky-20251001",
  "RefineSvGenotypesWithSnvs.min_ld_r2": 0.1,
  "RefineSvGenotypesWithSnvs.min_sv_ac": 50,
  "RefineSvGenotypesWithSnvs.min_sv_af": 0.001,
  "RefineSvGenotypesWithSnvs.min_an": 2000,
  "RefineSvGenotypesWithSnvs.output_prefix": "dfci-g2c.v1.\$CONTIG",
  "RefineSvGenotypesWithSnvs.sample_covariates": "$MAIN_WORKSPACE_BUCKET/data/sv_regenotyping/dfci-g2c.v1.sv_imputation_covariates.tsv.gz",
  "RefineSvGenotypesWithSnvs.sample_group_labels": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.qc_ancestry.tsv",
  "RefineSvGenotypesWithSnvs.snv_exclusion_bed": "$MAIN_WORKSPACE_BUCKET/data/sv_regenotyping/dfci-g2c.v1.sv_regenotyping.snv_mask.bed.gz",
  "RefineSvGenotypesWithSnvs.snv_freq_scalar": 10,
  "RefineSvGenotypesWithSnvs.snv_vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/data/sv_regenotyping/dfci-g2c.v1.sv_regenotyping.snv_vcf_info.\$CONTIG.tsv",
  "RefineSvGenotypesWithSnvs.sv_vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/\$CONTIG/HardFilterPart2/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.identical.reclustered.posthoc_filtered.vcf.gz",
  "RefineSvGenotypesWithSnvs.sv_vcf_idx": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/\$CONTIG/HardFilterPart2/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.identical.reclustered.posthoc_filtered.vcf.gz.tbi"
}
EOF

# Submit, monitor, and stage/cleanup SV GT refinement
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/RefineSvGenotypesWithSnvs.wdl \
  --input-json-template $staging_dir/RefineSvGenotypesWithSnvs.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup \
  --name RefineSvGenotypesWithSnvs \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.RefineSvGenotypesWithSnvs.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --max-attempts 3

# (Dev) chr19 submission
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/code/wdl/pancan_germline_wgs/RefineSvGenotypesWithSnvs.wdl \
  code/wdl/pancan_germline_wgs/ && \
cromshell --no_turtle -t 120 -mc submit \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip g2c.dependencies.zip \
  --no-validation code/wdl/pancan_germline_wgs/RefineSvGenotypesWithSnvs.wdl \
  /home/jupyter/cromshell/inputs/RefineSvGenotypesWithSnvs.inputs.chr19.json \
| jq .id | tr -d '"' > scratch/dev.wid && \
monitor_workflow $( cat scratch/dev.wid ) 1

# (Dev) debug
cromshell -t 180 metadata $( cat scratch/dev.wid ) | fgrep -A30 -B30 ailed | head -n500


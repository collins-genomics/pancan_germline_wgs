#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Initial quality control of G2C germline callset after joint genotyping

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


###########################################
# Compile .fam file of reported relatives #
###########################################

# Initialize staging directory
staging_dir=staging/fam_curation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Download sample metadata after GATK-HC
gsutil cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/

# Localize external .fam files staged in G2C bucket
gsutil cp -r \
  gs://dfci-g2c-refs/fam \
  $staging_dir/

# Curate 1000 Genomes project trios
code/scripts/external_fam_to_g2c.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --in-fam $staging_dir/fam/hgsvc.ped \
  --cohort hgsvc \
  --out-fam $staging_dir/dfci-g2c.hgsvc.fam

# Curate CEPH pedigrees
code/scripts/external_fam_to_g2c.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --in-fam $staging_dir/fam/ceph.remapped_ids.fam \
  --cohort ceph \
  --out-fam $staging_dir/dfci-g2c.ceph.fam

# Curate UFC pedigrees
code/scripts/external_fam_to_g2c.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --in-fam $staging_dir/fam/dfci_ufc_cases.fam \
  --cohort ufc \
  --out-fam $staging_dir/dfci-g2c.ufc.fam

# Curate MESA pedigrees
code/scripts/external_fam_to_g2c.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --in-fam $staging_dir/fam/mesa.fam \
  --cohort mesa \
  --out-fam $staging_dir/dfci-g2c.mesa.fam

# Collapse all cohort-specific .fam files
cat $staging_dir/dfci-g2c.*.fam \
| sort -Vk1,1 -k2,2V \
> $staging_dir/dfci-g2c.reported_families.fam

# Stage combined .fam file in Google bucket for future use
gsutil cp \
  $staging_dir/dfci-g2c.reported_families.fam \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/


###################################################
# Identify candidate twins / technical replicates #
###################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/fam_curation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Build list of all autosomal GATK-HC VCFs
for suf in vcfs tbis; do
  flist=$staging_dir/gatkhc.$suf.uris.list
  if [ -e $flist ]; then rm $flist; fi
done
for k in $( seq 1 22 ); do
  contig="chr$k"

  json_fname="PosthocCleanupPart2.$contig.outputs.json"
  gsutil -m cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart2/$contig/$json_fname \
    $staging_dir/
  
  jq .filtered_vcfs $staging_dir/$json_fname \
  | sed 's/,/\n/g' | tr -d '[]"' | sed '/^$/d' | awk '{ print $1 }' \
  >> $staging_dir/gatkhc.vcfs.uris.list
  
  jq .filtered_vcf_idxs $staging_dir/$json_fname \
  | sed 's/,/\n/g' | tr -d '[]"' | sed '/^$/d' | awk '{ print $1 }' \
  >> $staging_dir/gatkhc.tbis.uris.list
done

# Write input .json
cat << EOF | python -m json.tool > cromshell/inputs/InferTwins.inputs.json
{
  "InferTwins.ConcatVcfs.bcftools_concat_options": "--allow-overlaps",
  "InferTwins.ConcatVcfs.cpu_cores": 4,
  "InferTwins.ConcatVcfs.disk_gb": 1000,
  "InferTwins.ConcatVcfs.mem_gb": 16,
  "InferTwins.PrepVcf.mem_gb": 8.0,
  "InferTwins.PrepVcf.n_cpu": 4,
  "InferTwins.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "InferTwins.min_kin": 0.45,
  "InferTwins.only_snvs": true,
  "InferTwins.out_prefix": "dfci-g2c.v1",
  "InferTwins.plink2_docker": "vanallenlab/g2c_pipeline:plink2",
  "InferTwins.vcfs": $( collapse_txt $staging_dir/gatkhc.vcfs.uris.list ),
  "InferTwins.vcf_idxs": $( collapse_txt $staging_dir/gatkhc.tbis.uris.list )
}
EOF

# Submit twin inference workflow
cromshell --no_turtle -t 120 -mc submit --no-validation \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip g2c.dependencies.zip \
  code/wdl/pancan_germline_wgs/InferTwins.wdl \
  cromshell/inputs/InferTwins.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list )

# Once workflow is complete, stage output
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell-*/InferTwins/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage

# Locally filter putative twins to remove pairs that are more likely to be 
# duplicated WGS files (based on WGS QC metric similarity) than true replicates
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.kin0.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/
code/scripts/clean_candidate_twins.R \
  $staging_dir/dfci-g2c.v1.kin0.gz \
  $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/dfci-g2c.v1.cleaned.kin0
gzip -f $staging_dir/dfci-g2c.v1.cleaned.kin0
zcat $staging_dir/dfci-g2c.v1.cleaned.kin0.gz \
| fgrep -v "#" | cut -f2,4 \
> $staging_dir/dfci-g2c.v1.cleaned.tsv
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.cleaned.kin0.gz \
  $staging_dir/dfci-g2c.v1.cleaned.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/


#########################################################################
# Curate existing long-read and short-read WGS callsets for AoU samples #
#########################################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/aou_benchmark_curation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Collate list of AoU sample IDs corresponding to those samples present after variant calling
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| gunzip -c \
| awk -v FS="\t" -v OFS="\t" '{ if ($NF=="True" && $3=="aou") print $2 }' \
> $staging_dir/aou_ids.present_after_calling.samples.list
gsutil -m cp \
  $staging_dir/aou_ids.present_after_calling.samples.list \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/

# Exclude samples that failed srWGS snv/indel QC from srWGS consideration
gsutil -u $GPROJECT -m cat \
  gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv \
| sed '1d' | cut -f1 \
| fgrep -xvf - $staging_dir/aou_ids.present_after_calling.samples.list \
> $staging_dir/aou_ids.present_after_calling.no_srwgs_flagged.samples.list
gsutil -m cp \
  $staging_dir/aou_ids.present_after_calling.no_srwgs_flagged.samples.list \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/

# Find intersection of samples that have both short variants and SVs from lrWGS
# and were also present at the end of G2C variant calling
gsutil -u $GPROJECT -m cat \
  gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_vcf/GRCh38/cohort_for_GLNexus_2023Q1_1027.g.vcf.bgz \
| gunzip -c | head -n5000 | bcftools query -l \
> $staging_dir/aou.lrwgs.snv_indel.samples.list
gsutil -u $GPROJECT -m cat \
  gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_sv/GRCh38/integrated_sv_with_hprc_year_1_more_stringent.vcf.gz \
| gunzip -c | head -n5000 | bcftools query -l \
> $staging_dir/aou.lrwgs.sv.samples.list
fgrep -xf \
  $staging_dir/aou.lrwgs.snv_indel.samples.list \
  $staging_dir/aou.lrwgs.sv.samples.list \
| sort -V | uniq \
> $staging_dir/aou.lrwgs.has_complete_variation.samples.list
gsutil -m cp \
  $staging_dir/aou.lrwgs.has_complete_variation.samples.list \
  $MAIN_WORKSPACE_BUCKET/refs/aou/
fgrep -xf \
  $staging_dir/aou.lrwgs.has_complete_variation.samples.list \
  $staging_dir/aou_ids.present_after_calling.samples.list \
> $staging_dir/aou_ids.present_after_calling.complete_lrwgs.samples.list
gsutil -m cp \
  $staging_dir/aou_ids.present_after_calling.complete_lrwgs.samples.list \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/

# Write input .json for SV curation workflow
cat << EOF | python -m json.tool > cromshell/inputs/PreprocessAouSvs.inputs.json
{
 "PreprocessAouSvs.g2c_pipeline_docker": "vanallenlab/g2c_pipeline:050a724",
 "PreprocessAouSvs.srwgs_samples_list": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/aou_ids.present_after_calling.no_srwgs_flagged.samples.list",
 "PreprocessAouSvs.lrwgs_samples_list": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/aou_ids.present_after_calling.complete_lrwgs.samples.list"
}
EOF

# Submit SV curation workflow
cromshell --no_turtle -t 120 -mc submit --no-validation \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/external_data_curation/PreprocessAouSvs.wdl \
  cromshell/inputs/PreprocessAouSvs.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.PreprocessAouSvs.job_ids.list

# Monitor SV curation workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PreprocessAouSvs.job_ids.list ) 2

# Once workflow is complete, stage output
cromshell -t 120 list-outputs --json-summary \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PreprocessAouSvs.job_ids.list ) \
| python -m json.tool \
> $staging_dir/PreprocessAouSvs.outputs.json
cat \
  <( jq .\"PreprocessAouSvs.cleaned_lrwgs_vcf_idxs\" $staging_dir/PreprocessAouSvs.outputs.json ) \
  <( jq .\"PreprocessAouSvs.cleaned_lrwgs_vcfs\" $staging_dir/PreprocessAouSvs.outputs.json ) \
| awk '{ print $1 }' | tr -d '[]",' | sed '/^$/d' \
| gsutil -m cp -I $MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/
cat \
  <( jq .\"PreprocessAouSvs.cleaned_srwgs_vcf_idxs\" $staging_dir/PreprocessAouSvs.outputs.json ) \
  <( jq .\"PreprocessAouSvs.cleaned_srwgs_vcfs\" $staging_dir/PreprocessAouSvs.outputs.json ) \
| awk '{ print $1 }' | tr -d '[]",' | sed '/^$/d' \
| gsutil -m cp -I $MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/

# Clear Cromwell execution & output buckets for SV curation workflow
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.PreprocessAouSvs.job_ids.list \
                | awk \
                  -v exec_prefix="$WORKSPACE_BUCKET/cromwell-execution/PreprocessAouSvs/" \
                  -v out_prefix="$WORKSPACE_BUCKET/cromwell/outputs/PreprocessAouSvs/" \
                  -v OFS="\n" \
                  '{ print exec_prefix$1"/**", out_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage

# Write input .json for lrWGS SNV/indel curation
cat << EOF | python -m json.tool > cromshell/inputs/PreprocessAouLrwgsSnvs.inputs.json
{
 "PreprocessAouLrwgsSnvs.g2c_pipeline_docker": "vanallenlab/g2c_pipeline:7d94d38",
 "PreprocessAouLrwgsSnvs.ref_fasta" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
 "PreprocessAouLrwgsSnvs.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
 "PreprocessAouLrwgsSnvs.samples_list": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/aou_ids.present_after_calling.complete_lrwgs.samples.list"
}
EOF

# Submit lrWGS SNV/indel curation workflow
cromshell --no_turtle -t 120 -mc submit --no-validation \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/external_data_curation/PreprocessAouLrwgsSnvs.wdl \
  cromshell/inputs/PreprocessAouLrwgsSnvs.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.PreprocessAouLrwgsSnvs.job_ids.list

# Monitor lrWGS SNV/indel curation workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PreprocessAouLrwgsSnvs.job_ids.list ) 5

# Stage lrWGS short variants once curated
cromshell -t 120 list-outputs \
$( tail -n1 cromshell/job_ids/dfci-g2c.v1.PreprocessAouLrwgsSnvs.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/snv_indel/

# Clean up execution bucket for lrWGS short variant curation
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.PreprocessAouLrwgsSnvs.job_ids.list \
                | awk \
                  -v exec_prefix="$WORKSPACE_BUCKET/cromwell-execution/PreprocessAouLrwgsSnvs/" \
                  -v out_prefix="$WORKSPACE_BUCKET/cromwell/outputs/PreprocessAouLrwgsSnvs/" \
                  -v OFS="\n" \
                  '{ print exec_prefix$1"/**", out_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage

# Gather list of AoU samples that had SV data (this is not the entire cohort)
gsutil -m cat \
  gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.chr22.vcf.gz \
| gunzip -c | head -n5000 | bcftools query -l \
> $staging_dir/AoU.G2C_samples_with_srWGS_sv.aou_ids.list
gsutil -m cp \
  $staging_dir/AoU.G2C_samples_with_srWGS_sv.aou_ids.list \
  $MAIN_WORKSPACE_BUCKET/refs/aou/

# Write desired header for srWGS SNV/indel output VCFs
gsutil -m cat \
  gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/snv_indel/1KGP.srWGS.snv_indel.cleaned.chr1.vcf.gz \
| gunzip -c | head -n5000 | fgrep "##" \
> $staging_dir/AoU.srwgs.snv_indel_header.vcf

# Important note: raw AoU srWGS SNV/indel curation is documented elsewhere
# This is due to the need to manipulate the raw callset in Hail .vds format,
# which needs to be run on a dedicated Hail dataproc spark cluster
# See: pancan_germline_wgs/shell/aou_rw/extract_aou_srwgs_short_variants.sh

# Index AoU srWGS short variants after extraction from the overall VDS (see above)

# Write template input .json for srWGS short variant indexing
cat << EOF > $staging_dir/IndexAouSrwgsSnvs.inputs.template.json
{
  "IndexVcf.copy_index_to_vcf_bucket": true,
  "IndexVcf.vcf": "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/snv_indel/AoU.srWGS.snv_indel.cleaned.\$CONTIG.vcf.bgz"
}
EOF

# Submit, monitor, stage, and cleanup AoU srWGS short variant indexing workflows
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/IndexVcf.wdl \
  --input-json-template $staging_dir/IndexAouSrwgsSnvs.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/scratch/ \
  --name IndexAouSrwgsSnvs \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.IndexAouSrwgsSnvs.initial_qc.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --max-attempts 4


###############################################################
# Build QC sample inclusion priority & sampling probabilities #
###############################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/sample_priority
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Localize most up-to-date sample QC manifest
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/

# Generate sampling probabilities
code/scripts/assign_sample_qc_weights.R \
  --qc-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --pass-column filtersites_qc_pass \
  --pass-column gatksv_posthoc_qc_pass \
  --pass-column gatkhc_posthoc_qc_pass \
  --out-tsv $staging_dir/dfci-g2c.v1.qc.sample_weights.tsv

# Extract list of sample IDs still present in callset at this stage
zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| awk -v FS="\t" -v OFS="\t" '{ if ($NF=="True") print $1 }' \
> $staging_dir/g2c_ids.present_after_calling.samples.list
zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| awk -v FS="\t" -v OFS="\t" '{ if ($NF=="True") print $1, $2, $3 }' \
> $staging_dir/g2c_ids.present_after_calling.sample_cohort_map.tsv

# Find list of samples included in complete trios
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam \
  $staging_dir/
code/scripts/subset_fam.R \
  --in-fam $staging_dir/dfci-g2c.reported_families.fam \
  --all-samples $staging_dir/g2c_ids.present_after_calling.samples.list \
  --out-fam $staging_dir/dfci-g2c.reported_families.filtered.fam
awk -v FS="\t" -v OFS="\t" '{ if ($2!=0 && $3!=0 && $4!=0) print }' \
  $staging_dir/dfci-g2c.reported_families.filtered.fam \
> $staging_dir/dfci-g2c.reported_families.filtered.complete.fam
awk -v OFS="\n" '{ print $2, $3, $4 }' \
  $staging_dir/dfci-g2c.reported_families.filtered.complete.fam \
| sort -V | uniq \
> $staging_dir/dfci-g2c.reported_families.filtered.complete.samples.list

# Find list of samples included in complete twin pairs
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.kin0.gz \
  $staging_dir/
zcat $staging_dir/dfci-g2c.v1.cleaned.kin0.gz \
| fgrep -v "#" | awk -v OFS="\t" '{ print "twins_"NR, $2, $4, 0, 0 }' \
> $staging_dir/all_twins.fam
code/scripts/subset_fam.R \
  --in-fam $staging_dir/all_twins.fam \
  --all-samples $staging_dir/g2c_ids.present_after_calling.samples.list \
  --out-fam $staging_dir/all_twins.filtered.fam
awk -v FS="\t" -v OFS="\n" '{ if ($2!=0 && $3!=0) print $2, $3 }' \
  $staging_dir/all_twins.filtered.fam \
| sort -V | uniq \
> $staging_dir/dfci-g2c.inferred_twins.complete.samples.list

# Get list of samples with complete HGSV long-read WGS calls
# Supplement with twins/replicates, since these are equally useful for our purposes
gsutil -m cat \
  gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/snv_indel/1KGP.lrWGS.snv_indel.cleaned.chrY.vcf.gz \
| gunzip -c | head -n10000 | fgrep "#" | bcftools query -l \
> $staging_dir/1KGP.lrWGS.snv_indel.external_ids.tsv
gsutil -m cat \
  gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.chrY.vcf.gz \
| gunzip -c | head -n10000 | fgrep "#" | bcftools query -l \
> $staging_dir/1KGP.lrWGS.sv.external_ids.tsv
fgrep -xf \
  $staging_dir/1KGP.lrWGS.snv_indel.external_ids.tsv \
  $staging_dir/1KGP.lrWGS.sv.external_ids.tsv \
| sort | uniq \
> $staging_dir/1KGP.lrWGS.complete.external_ids.tsv
awk -v OFS="\t" '{ if ($3=="hgsvc") print $2, $1 }' \
  $staging_dir/g2c_ids.present_after_calling.sample_cohort_map.tsv \
| sort -k1,1 \
| join -j 1 -t $'\t' - $staging_dir/1KGP.lrWGS.complete.external_ids.tsv \
| cut -f2 | sort -V \
> $staging_dir/1KGP.lrWGS.complete.g2c_ids.no_twins.samples.list
zcat $staging_dir/dfci-g2c.v1.cleaned.kin0.gz | cut -f2,4 \
| fgrep -wf $staging_dir/1KGP.lrWGS.complete.g2c_ids.no_twins.samples.list \
| sed 's/\t/\n/g' \
| cat - $staging_dir/1KGP.lrWGS.complete.g2c_ids.no_twins.samples.list \
| sort -V | uniq \
> $staging_dir/1KGP.lrWGS.complete.g2c_ids.samples.list

# Get list of samples with complete AoU long-read WGS calls
# Supplement with twins/replicates, since these are equally useful for our purposes
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/snv_indel/AoU.lrWGS.snv_indel.cleaned.chrY.vcf.gz \
| gunzip -c | head -n10000 | fgrep "#" | bcftools query -l \
> $staging_dir/AoU.lrWGS.snv_indel.external_ids.tsv
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.chrY.vcf.gz \
| gunzip -c | head -n10000 | fgrep "#" | bcftools query -l \
> $staging_dir/AoU.lrWGS.sv.external_ids.tsv
fgrep -xf \
  $staging_dir/AoU.lrWGS.snv_indel.external_ids.tsv \
  $staging_dir/AoU.lrWGS.sv.external_ids.tsv \
| sort | uniq \
> $staging_dir/AoU.lrWGS.complete.external_ids.tsv
awk -v OFS="\t" '{ if ($3=="aou") print $2, $1 }' \
  $staging_dir/g2c_ids.present_after_calling.sample_cohort_map.tsv \
| sort -k1,1 \
| join -j 1 -t $'\t' - $staging_dir/AoU.lrWGS.complete.external_ids.tsv \
| cut -f2 | sort -V \
> $staging_dir/AoU.lrWGS.complete.g2c_ids.no_twins.samples.list
zcat $staging_dir/dfci-g2c.v1.cleaned.kin0.gz | cut -f2,4 \
| fgrep -wf $staging_dir/AoU.lrWGS.complete.g2c_ids.no_twins.samples.list \
| sed 's/\t/\n/g' \
| cat - $staging_dir/AoU.lrWGS.complete.g2c_ids.no_twins.samples.list \
| sort -V | uniq \
> $staging_dir/AoU.lrWGS.complete.g2c_ids.samples.list

# Combine lists of samples with lrWGS available from at least one source
cat \
  $staging_dir/1KGP.lrWGS.complete.g2c_ids.samples.list \
  $staging_dir/AoU.lrWGS.complete.g2c_ids.samples.list \
| sort -V | uniq \
> $staging_dir/dfci-g2c.lrWGS.complete.samples.list

# Compute sample priority and append sampling probabilities
cat \
  $staging_dir/dfci-g2c.reported_families.filtered.complete.samples.list \
  $staging_dir/dfci-g2c.inferred_twins.complete.samples.list \
  $staging_dir/dfci-g2c.lrWGS.complete.samples.list \
  $staging_dir/g2c_ids.present_after_calling.samples.list \
| sort -V | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | sort -k1,1 \
| join -j 1 -t $'\t' \
  - <( sort -k1,1 $staging_dir/dfci-g2c.v1.qc.sample_weights.tsv ) \
| sort -nrk2,2 -k3,3nr -k1,1V \
| cat <( echo -e "#G2C_ID\tpriority\tweight" ) - \
> $staging_dir/dfci-g2c.v1.sample_qc_priority.tsv
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.sample_qc_priority.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/


#####################################
# Build inter-cohort sample ID maps #
#####################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/external_id_maps
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Localize most up-to-date sample QC manifest & cleaned twins
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.kin0.gz \
  $staging_dir/

# Build map of sample IDs for 1KGP (accounting for external twins/replicates)
code/scripts/make_id_map.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --cohort hgsvc \
  --twins $staging_dir/dfci-g2c.v1.cleaned.kin0.gz \
  --out-tsv $staging_dir/dfci-g2c.v1.1KGP_id_map.tsv

# Build map of sample IDs for AoU (accounting for external twins/replicates)
code/scripts/make_id_map.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --cohort aou \
  --twins $staging_dir/dfci-g2c.v1.cleaned.kin0.gz \
  --out-tsv $staging_dir/dfci-g2c.v1.AoU_id_map.tsv

# Copy ID maps to reference directory
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.1KGP_id_map.tsv \
  $staging_dir/dfci-g2c.v1.AoU_id_map.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/


########################################################
# Build ancestry & case:control labels for QC plotting #
########################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/ancestry_phenotype_maps
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Localize most up-to-date sample QC manifest & cleaned twins
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/

# Extract ancestry map
cidx=$( zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
        | head -n1 | sed 's/\t/\n/g' \
        | awk -v FS="\t" '{ if ($1=="intake_qc_pop") print NR }' )
zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| cut -f1,$cidx | grep -ve '^G2C_id' | sort -Vk1,1 \
> $staging_dir/dfci-g2c.v1.qc_ancestry.tsv

# Extract case|control map, reassigning unknowns to controls
cidx=$( zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
        | head -n1 | sed 's/\t/\n/g' \
        | awk -v FS="\t" '{ if ($1=="batching_pheno") print NR }' )
zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| cut -f1,$cidx | grep -ve '^G2C_id' | sort -Vk1,1 \
> $staging_dir/dfci-g2c.v1.qc_phenotype.tsv

# Copy ID maps to reference directory
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.qc_ancestry.tsv \
  $staging_dir/dfci-g2c.v1.qc_phenotype.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/


##############################
# Collect initial QC metrics #
##############################

# Note: this workflow below is scattered across all five workspaces for 
# max parallelization. It must be submitted as below in each workspace.

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Check to ensure there is a local copy of calling intervals
if ! [ -e $staging_dir/calling_intervals ]; then
  mkdir $staging_dir/calling_intervals
  gsutil -m cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/*.sharded.interval_list \
    $staging_dir/calling_intervals/
fi

# Write two-column .tsv of VCF & index info for each contig
while read contig; do
  gsutil cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart2/$contig/PosthocCleanupPart2.$contig.outputs.json \
  | jq '.["PosthocCleanupPart2.filtered_vcfs"]' \
  | fgrep "gs://" | awk '{ print $1 }' | tr -d '",' \
  | cat - <( echo -e "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/$contig/HardFilterPart2/dfci-g2c.v1.$contig.concordance.gq_recalibrated.identical.reclustered.posthoc_filtered.vcf.gz" ) \
  | awk -v OFS="\t" '{ print $1, $1".tbi" }' \
  > $staging_dir/dfci-g2c.v1.initial_qc.vcf_info.$contig.tsv
done < contig_lists/dfci-g2c.v1.contigs.$WN.list
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.initial_qc.vcf_info.*.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/vcf_list_inputs/

# Initialize .json of contig-specific overrides for scatter counts
echo "{ " > $staging_dir/CollectInitialVcfQcMetrics.contig_variable_overrides.json
while read contig; do
  kc=$( fgrep -v "@" \
          $staging_dir/calling_intervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list \
        | wc -l | awk '{ printf "%i\n", $1 / 3 }' )
  echo "\"$contig\" : {\"CONTIG_SCATTER_COUNT\" : $kc},"
done < contig_lists/dfci-g2c.v1.contigs.$WN.list \
| paste -s -d\  | sed 's/,$//g' \
>> $staging_dir/CollectInitialVcfQcMetrics.contig_variable_overrides.json
echo " }" >> $staging_dir/CollectInitialVcfQcMetrics.contig_variable_overrides.json

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectInitialVcfQcMetrics.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": \$CONTIG_SCATTER_COUNT,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.BenchmarkSites.indel_mem_scalar": 2.0,
  "CollectVcfQcMetrics.BenchmarkSites.snv_mem_scalar": 4.0,
  "CollectVcfQcMetrics.ChunkCommonVcf.mem_gb": 7.5,
  "CollectVcfQcMetrics.ChunkCommonVcf.cpu_cores": 4,
  "CollectVcfQcMetrics.common_af_cutoff": 0.001,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:1aac84d",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 5000,
  "CollectVcfQcMetrics.output_prefix": "dfci-g2c.v1.initial_qc.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/snv_indel/1KGP.srWGS.snv_indel.cleaned.\$CONTIG.vcf.gz",
                                                 "gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/snv_indel/AoU.srWGS.snv_indel.cleaned.\$CONTIG.vcf.bgz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/snv_indel/1KGP.lrWGS.snv_indel.cleaned.\$CONTIG.vcf.gz",
                                                 "gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/snv_indel/AoU.lrWGS.snv_indel.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/snv_indel/1KGP.srWGS.snv_indel.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/snv_indel/AoU.srWGS.snv_indel.cleaned.\$CONTIG.vcf.bgz",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/snv_indel/1KGP.lrWGS.snv_indel.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/snv_indel/AoU.lrWGS.snv_indel.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.sample_priority_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.sample_qc_priority.tsv",
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.snv.sites.bed.gz"],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.indel.sites.bed.gz"],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/vcf_list_inputs/dfci-g2c.v1.initial_qc.vcf_info.\$CONTIG.tsv"
}
EOF

# Submit, monitor, stage, and cleanup QC metric collection workflows
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectInitialVcfQcMetrics.inputs.template.json \
  --contig-variable-overrides $staging_dir/CollectInitialVcfQcMetrics.contig_variable_overrides.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/VcfQcMetrics/ \
  --name CollectInitialVcfQcMetrics \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.CollectInitialVcfQcMetrics.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 60 \
  --vm-gate 400 \
  --submission-gate 60 \
  --max-attempts 3


#####################
# Curate QC targets #
#####################

# Reaffirm staging directory
staging_dir=staging/qc_targets
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Estimate number of variants per genome in gnomAD for the necessary contigs
for k in $( seq 1 22 ) X Y; do
  gsutil cat \
    gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/chr$k/gnomad.v4.1.chr$k.*.sites.bed.gz
done | gunzip -c \
| code/scripts/estimate_vpg_from_sites.py \
| fgrep -v "#" \
| awk -v OFS="\t" '{ print "variants_per_genome."$1":median", $2 }' \
> $staging_dir/dfci-g2c.v1.qc_targets.tsv

# Copy QC targets to central bucket for reference by Cromwell
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.qc_targets.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/


##########################################
# Analyze & visualize initial QC metrics #
##########################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

cat << EOF > $staging_dir/main_keys.list
size_distrib
af_distrib
size_vs_af_distrib
all_svs_bed
common_snvs_bed
common_indels_bed
common_svs_bed
genotype_distrib
ld_stats
EOF

cat << EOF > $staging_dir/bench_keys.list
site_benchmark_ppv_by_freqs
site_benchmark_sensitivity_by_freqs
site_benchmark_common_snv_ppv_beds
site_benchmark_common_indel_ppv_beds
site_benchmark_common_sv_ppv_beds
site_benchmark_common_snv_sens_beds
site_benchmark_common_indel_sens_beds
site_benchmark_common_sv_sens_beds
twin_genotype_benchmark_distribs
trio_mendelian_violation_distribs
EOF

# Clear old input arrays
while read key; do  
  fname=$staging_dir/$key.uris.list
  if [ -e $fname ]; then rm $fname; fi
done < $staging_dir/main_keys.list
while read key; do
  for subset in giab_easy giab_hard; do
    fname=$staging_dir/$key.$subset.uris.list
    if [ -e $fname ]; then rm $fname; fi
  done
done < $staging_dir/bench_keys.list
for suffix in af_distribution size_distribution; do
  fname=$staging_dir/gnomAD_$suffix.uris.list
  if [ -e $fname ]; then rm $fname; fi
done
for key in sample_benchmark_ppv_distribs sample_benchmark_sensitivity_distribs; do
  for dset in external_srwgs external_lrwgs; do
    for subset in giab_easy giab_hard; do
      fname=$staging_dir/$key.$subset.$dset.uris.list
      if [ -e $fname ]; then rm $fname; fi
    done
  done
done

# Build input arrays
for k in $( seq 1 22 ) X Y; do
  
  # Localize output tracker json and get URIs for QC metrics
  json_fname=CollectInitialVcfQcMetrics.chr$k.outputs.json
  gsutil cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/VcfQcMetrics/chr$k/$json_fname \
    $staging_dir/
  while read key; do
    jq .\"CollectVcfQcMetrics.$key\" $staging_dir/$json_fname \
    | fgrep -xv "null" | tr -d '"' \
    >> $staging_dir/$key.uris.list
  done < $staging_dir/main_keys.list
  while read key; do
    for subset in giab_easy giab_hard; do
      jq .\"CollectVcfQcMetrics.$key\" $staging_dir/$json_fname \
      | fgrep -xv "null" | tr -d '"[]' | sed 's/,$/\n/g' \
      | sed '/^$/d' | awk '{ print $1 }' | fgrep $subset \
      >> $staging_dir/$key.$subset.uris.list
    done
  done < $staging_dir/bench_keys.list

  # Due to delisting behavior of manage_chromshards.py, external sample benchmark
  # results need to be parsed in a custom manner as below
  for key in sample_benchmark_ppv_distribs sample_benchmark_sensitivity_distribs; do
    for dset in external_srwgs external_lrwgs; do
      for subset in giab_easy giab_hard; do
        jq .\"CollectVcfQcMetrics.$key\" $staging_dir/$json_fname \
        | fgrep -xv "null" | tr -d '"[]' | sed 's/,$/\n/g' \
        | sed '/^$/d' | awk '{ print $1 }' | fgrep $dset | fgrep $subset \
        >> $staging_dir/$key.$subset.$dset.uris.list
      done
    done
  done

  # Clear local copy of output tracker json
  rm $staging_dir/$json_fname

  # Add precomputed gnomAD v4.1 reference distributions to file lists
  for suffix in af_distribution size_distribution; do
    echo "gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/chr$k/gnomad.v4.1.chr$k.$suffix.merged.tsv.gz" \
    >> $staging_dir/gnomAD_$suffix.uris.list
  done
done

# Write input .json
cat << EOF | python -m json.tool > cromshell/inputs/PlotInitialVcfQcMetrics.inputs.json
{
  "PlotVcfQcMetrics.af_distribution_tsvs": $( collapse_txt $staging_dir/af_distrib.uris.list ),
  "PlotVcfQcMetrics.all_sv_beds": $( collapse_txt $staging_dir/all_svs_bed.uris.list ),
  "PlotVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PlotVcfQcMetrics.benchmark_interval_names": ["Easy", "Hard"],
  "PlotVcfQcMetrics.common_af_cutoff": 0.001,
  "PlotVcfQcMetrics.common_snv_beds": $( collapse_txt $staging_dir/common_snvs_bed.uris.list ),
  "PlotVcfQcMetrics.common_indel_beds": $( collapse_txt $staging_dir/common_indels_bed.uris.list ),
  "PlotVcfQcMetrics.common_sv_beds": $( collapse_txt $staging_dir/common_svs_bed.uris.list ),
  "PlotVcfQcMetrics.custom_qc_target_metrics": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.qc_targets.tsv",
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:1aac84d",
  "PlotVcfQcMetrics.output_prefix": "dfci-g2c.v1.initial_qc",
  "PlotVcfQcMetrics.peak_ld_stat_tsvs": $( collapse_txt $staging_dir/ld_stats.uris.list ),
  "PlotVcfQcMetrics.PlotSiteBenchmarking.mem_gb": 32,
  "PlotVcfQcMetrics.PlotSiteBenchmarking.n_cpu": 8,
  "PlotVcfQcMetrics.PlotSiteMetrics.mem_gb": 32,
  "PlotVcfQcMetrics.PlotSiteMetrics.n_cpu": 8,
  "PlotVcfQcMetrics.ref_af_distribution_tsvs": $( collapse_txt $staging_dir/gnomAD_af_distribution.uris.list ),
  "PlotVcfQcMetrics.ref_size_distribution_tsvs": $( collapse_txt $staging_dir/gnomAD_size_distribution.uris.list ),
  "PlotVcfQcMetrics.ref_cohort_prefix": "gnomAD_v4.1",
  "PlotVcfQcMetrics.ref_cohort_plot_title": "gnomAD v4.1",
  "PlotVcfQcMetrics.sample_ancestry_labels": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.qc_ancestry.tsv",
  "PlotVcfQcMetrics.sample_phenotype_labels": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.qc_phenotype.tsv",
  "PlotVcfQcMetrics.sample_benchmark_dataset_prefixes": ["external_srwgs", "external_lrwgs"],
  "PlotVcfQcMetrics.sample_benchmark_dataset_titles": ["external srWGS", "external lrWGS"],
  "PlotVcfQcMetrics.sample_benchmark_ppv_distribs": [[ $( collapse_txt $staging_dir/sample_benchmark_ppv_distribs.giab_easy.external_srwgs.uris.list ),
                                                       $( collapse_txt $staging_dir/sample_benchmark_ppv_distribs.giab_hard.external_srwgs.uris.list ) ],
                                                     [ $( collapse_txt $staging_dir/sample_benchmark_ppv_distribs.giab_easy.external_lrwgs.uris.list ),
                                                       $( collapse_txt $staging_dir/sample_benchmark_ppv_distribs.giab_hard.external_lrwgs.uris.list ) ]],
  "PlotVcfQcMetrics.sample_benchmark_sensitivity_distribs": [[ $( collapse_txt $staging_dir/sample_benchmark_sensitivity_distribs.giab_easy.external_srwgs.uris.list ),
                                                               $( collapse_txt $staging_dir/sample_benchmark_sensitivity_distribs.giab_hard.external_srwgs.uris.list ) ],
                                                             [ $( collapse_txt $staging_dir/sample_benchmark_sensitivity_distribs.giab_easy.external_lrwgs.uris.list ),
                                                               $( collapse_txt $staging_dir/sample_benchmark_sensitivity_distribs.giab_hard.external_lrwgs.uris.list ) ]],
  "PlotVcfQcMetrics.sample_genotype_distribution_tsvs": $( collapse_txt $staging_dir/genotype_distrib.uris.list ),
  "PlotVcfQcMetrics.site_benchmark_common_snv_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_snv_ppv_beds.giab_easy.uris.list ),
                                                            $( collapse_txt $staging_dir/site_benchmark_common_snv_ppv_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_common_indel_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_indel_ppv_beds.giab_easy.uris.list ),
                                                              $( collapse_txt $staging_dir/site_benchmark_common_indel_ppv_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_common_sv_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_easy.uris.list ),
                                                           $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_common_snv_sens_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_snv_sens_beds.giab_easy.uris.list ),
                                                             $( collapse_txt $staging_dir/site_benchmark_common_snv_sens_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_common_indel_sens_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_indel_sens_beds.giab_easy.uris.list ),
                                                               $( collapse_txt $staging_dir/site_benchmark_common_indel_sens_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_common_sv_sens_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_sv_sens_beds.giab_easy.uris.list ),
                                                            $( collapse_txt $staging_dir/site_benchmark_common_sv_sens_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_ppv_by_freqs": [[ $( collapse_txt $staging_dir/site_benchmark_ppv_by_freqs.giab_easy.uris.list ),
                                                     $( collapse_txt $staging_dir/site_benchmark_ppv_by_freqs.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_sensitivity_by_freqs": [[ $( collapse_txt $staging_dir/site_benchmark_sensitivity_by_freqs.giab_easy.uris.list ),
                                                             $( collapse_txt $staging_dir/site_benchmark_sensitivity_by_freqs.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_dataset_prefixes": ["gnomad_v4.1"],
  "PlotVcfQcMetrics.site_benchmark_dataset_titles": ["gnomAD v4.1"],
  "PlotVcfQcMetrics.size_distribution_tsvs": $( collapse_txt $staging_dir/size_distrib.uris.list ),
  "PlotVcfQcMetrics.size_vs_af_distribution_tsvs": $( collapse_txt $staging_dir/size_vs_af_distrib.uris.list ),
  "PlotVcfQcMetrics.trio_mendelian_violation_distribs": [ $( collapse_txt $staging_dir/trio_mendelian_violation_distribs.giab_easy.uris.list ),
                                                          $( collapse_txt $staging_dir/trio_mendelian_violation_distribs.giab_hard.uris.list ) ],
  "PlotVcfQcMetrics.twin_genotype_benchmark_distribs": [ $( collapse_txt $staging_dir/twin_genotype_benchmark_distribs.giab_easy.uris.list ),
                                                         $( collapse_txt $staging_dir/twin_genotype_benchmark_distribs.giab_hard.uris.list ) ]
}
EOF

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotInitialVcfQcMetrics.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.PlotInitialVcfQcMetrics.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotInitialVcfQcMetrics.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotInitialVcfQcMetrics.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotQc/

# Clear Cromwell execution & output buckets for plotting job
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.PlotInitialVcfQcMetrics.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell/*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


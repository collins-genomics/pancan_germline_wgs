#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-HC joint genotyping on G2C phase 1

# Note that this code is designed to be run inside the AoU Researcher Workbench


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in scratch staging cromshell cromshell/inputs cromshell/job_ids cromshell/progress; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}
find code/ -name "*.R" | xargs -I {} chmod a+x {}

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh


# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Create dependencies .zip for G2C workflow submissions
cd code/wdl/pancan_germline_wgs && \
zip g2c.dependencies.zip *.wdl && \
mv g2c.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for GATK-HC workflow submissions
mkdir ~/scratch/gatkhc.dependencies && \
cp code/wdl/pancan_germline_wgs/Utilities.wdl ~/scratch/gatkhc.dependencies/ && \
cp code/wdl/gatk-hc/*.wdl ~/scratch/gatkhc.dependencies/ && \
cd ~/scratch/gatkhc.dependencies && \
zip gatkhc.dependencies.zip *.wdl && \
mv gatkhc.dependencies.zip ~/ && \
cd ~ && \
rm -rf ~/scratch/gatkhc.dependencies

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Ensure Cromwell/Cromshell are configured
code/scripts/setup_cromshell.py

# Install necessary packages
. code/refs/install_packages.sh python

# Download workspace-specific contig lists
gsutil cp -r \
  gs://dfci-g2c-refs/hg38/contig_lists \
  ./


#######################
# Generate sample map #
#######################

# Note: this only needs to be run once across all workspaces

# Refresh staging directory
staging_dir=staging/GenerateSampleMap
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Gather information needed for all samples to be included for joint genotyping
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/dfci-g2c.sample_meta.posthoc_outliers.tsv.gz \
| gunzip -c \
| awk -v FS="\t" -v OFS="\t" '{ if ($NF=="True") print $1, $2, $3 }' \
| sort -Vk1,1 \
> $staging_dir/dfci-g2c.v1.gatkhc.input_samples.info.tsv

# Write GATK-HC sample map
while read sid oid cohort; do
  if [ $cohort == "aou" ]; then
    uri_prefix=$MAIN_WORKSPACE_BUCKET
  else
    uri_prefix="gs:/"
  fi
  echo -e "$sid\t$uri_prefix/dfci-g2c-inputs/$cohort/gatk-hc/reblocked/$oid.reblocked.g.vcf.gz"
done < $staging_dir/dfci-g2c.v1.gatkhc.input_samples.info.tsv \
> $staging_dir/dfci-g2c.v1.gatkhc.sample_map.tsv

# Move sample map to main bucket for Cromwell access
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.gatkhc.sample_map.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/


#####################
# Prepare intervals #
#####################

# Note: this must be run once in each workspace

# Refresh staging directory
staging_dir=staging/PrepIntervals
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Download Broad standard hg38 calling intervals list
gsutil -m cp \
  gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list \
  $staging_dir/

# Split intervals into one per primary chromosome
ilist=$staging_dir/wgs_calling_regions.hg38.interval_list
for k in $( seq 1 22 ) X Y; do
  fgrep "@" $ilist > $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list
  awk -v contig="chr$k" '{ if ($1==contig) print }' $ilist \
  >> $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list
done

# Estimate ideal number of shards per chromosome
# Math as follows: 2,900 task quota x 5 workspaces = 14,500 total shards
# We will increase this value by 50% to ensure we aren't wasting much quota
# without also overloading the cromwell server
# Shards per chrom = 14,500 x 1.5 x chromosome size / genome size
genome_size=$( cat $staging_dir/gatkhc.wgs_calling_regions.hg38.chr*.interval_list \
               | fgrep -v "@" | cut -f1-3 | bedtools merge -i - \
               | awk '{ sum+=$3-$2 }END{ printf "%i\n", sum }' )

# Re-shard GATK-HC intervals per chromosome
while read contig; do
  # Get size of contig
  contig_size=$( cat $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.interval_list \
                 | fgrep -v "@" | cut -f1-3 | bedtools merge -i - \
                 | awk '{ sum+=$3-$2 }END{ printf "%i\n", sum }' )

  # Calculate number of desired shards
  n_shards=$( echo "14500" \
              | awk -v denom=$genome_size -v numer=$contig_size \
                '{ printf "%i\n", $1 * 1.5 * numer / denom }' )

  # Download gnomAD variant sites from this contig
  gsutil -m cp \
    gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/$contig/gnomad.v4.1.$contig.*.sites.bed.gz* \
    $staging_dir/
  find $staging_dir -name "*bed.gz" | xargs -I {} tabix -f {}

  # Estimate ideal number of variants per shard
  n_vars=$( zcat $staging_dir/gnomad.v4.1.$contig.*.sites.bed.gz | wc -l | awk '{ print $1-3 }' )
  vars_per_shard=$( echo $n_vars | awk -v denom=$n_shards '{ printf "%i\n", $1 / denom }' )

  # Shard GATK-HC intervals according to parameters determined above
  code/scripts/split_intervals.py \
    -i $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.interval_list \
    --var-sites $staging_dir/gnomad.v4.1.$contig.snv.sites.bed.gz \
    --var-sites $staging_dir/gnomad.v4.1.$contig.indel.sites.bed.gz \
    --var-sites $staging_dir/gnomad.v4.1.$contig.sv.sites.bed.gz \
    --vars-per-shard $vars_per_shard \
    --verbose \
    -o $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list

  # Clean up
  rm $staging_dir/gnomad.v4.1.*.*.sites.bed.gz*

done < contig_lists/dfci-g2c.v1.contigs.$WN.list

# Copy all sharded interval lists to google bucket for Cromwell access
gsutil cp \
  $staging_dir/gatkhc.wgs_calling_regions.hg38.chr*.sharded.interval_list \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/


####################
# Joint genotyping #
####################

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Refresh staging directory
staging_dir=staging/JointGenotyping
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Check to ensure there is a local copy of calling intervals
if ! [ -e staging/PrepIntervals ]; then
  mkdir staging/PrepIntervals
  gsutil -m cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/*.sharded.interval_list \
    staging/PrepIntervals/
fi

# Write .json of contig-specific scatter counts
echo "{ " > $staging_dir/contig_variable_overrides.json
while read contig; do
  kc=$( fgrep -v "@" \
          staging/PrepIntervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list \
        | wc -l | awk '{ printf "%i\n", $1 }' )
  echo "\"$contig\" : {\"CONTIG_SCATTER_COUNT\" : $kc },"
done < contig_lists/dfci-g2c.v1.contigs.$WN.list \
| paste -s -d\  | sed 's/,$//g' \
>> $staging_dir/contig_variable_overrides.json
echo " }" >> $staging_dir/contig_variable_overrides.json

# Write template .json for input
cat << EOF > $staging_dir/GnarlyJointGenotypingPart1.inputs.template.json
{
  "GnarlyJointGenotypingPart1.callset_name": "dfci-g2c.v1.\$CONTIG",
  "GnarlyJointGenotypingPart1.dbsnp_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
  "GnarlyJointGenotypingPart1.GnarlyGenotyperFT.machine_mem_mb": 12000,
  "GnarlyJointGenotypingPart1.gnarly_scatter_count": 1,
  "GnarlyJointGenotypingPart1.import_gvcfs_batch_size": 100,
  "GnarlyJointGenotypingPart1.import_gvcfs_disk_gb": 20,
  "GnarlyJointGenotypingPart1.ImportGVCFsFT.jvm_max_mb": 25000,
  "GnarlyJointGenotypingPart1.ImportGVCFsFT.machine_mem_mb": 44000,
  "GnarlyJointGenotypingPart1.ImportGVCFsFT.n_preemptible_tries": 0,
  "GnarlyJointGenotypingPart1.intervals_already_split": true,
  "GnarlyJointGenotypingPart1.make_hard_filtered_sites": false,
  "GnarlyJointGenotypingPart1.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
  "GnarlyJointGenotypingPart1.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "GnarlyJointGenotypingPart1.ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "GnarlyJointGenotypingPart1.sample_name_map": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/dfci-g2c.v1.gatkhc.sample_map.tsv",
  "GnarlyJointGenotypingPart1.top_level_scatter_count": \$CONTIG_SCATTER_COUNT,
  "GnarlyJointGenotypingPart1.unpadded_intervals_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/gatkhc.wgs_calling_regions.hg38.\$CONTIG.sharded.interval_list"
}
EOF

# Due to the multi-stage nature of joint genotyping (see below),
# we need to initialize the tracker .tsv where any contig with a bucket created
# in the parent staging bucket is considered to be staged.
# This is necessary because contigs can be effectively ~complete through this first
# phase of joint genotyping but persistently fail due to problematic loci/shards.
gsutil ls $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/ \
| xargs -I {} basename {} \
| awk -v OFS="\t" '{ print $1, "staged" }' \
> cromshell/progress/dfci-g2c.v1.JointGenotyping.progress.tsv

# Joint genotype per chromosome using chromsharded manager
# In practice, this workflow basically never completes successfully and needs
# to be handled in nested / iterative patches (see below)
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-hc/GnarlyJointGenotypingPart1.wdl \
  --input-json-template $staging_dir/GnarlyJointGenotypingPart1.inputs.template.json \
  --contig-variable-overrides $staging_dir/contig_variable_overrides.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/ \
  --dependencies-zip gatkhc.dependencies.zip \
  --name JointGenotyping \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.JointGenotyping.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 60 \
  --submission-gate 60 \
  --gcp-report-period 10 \
  --vm-gate 1000 \
  --max-attempts 2


################################################
# Semi-manual cleanup of genotyping first pass #
################################################

# Due to the piecemeal nature of how we ran joint genotyping in practice, we first
# need to collect all complete shards from the process above in a semi-manual
# manner before determining which intervals still need to be resharded/rerun

# This must be run once for each workspace

# Reaffirm staging directory
staging_dir=staging/JGSecondPass
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# First, copy all VCFs that completed during the first pass of joint genotyping
while read contig; do
  # Do nothing if contig was one of the development contigs processed earlier
  if [ $( gsutil ls $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart1/$contig 2>/dev/null | wc -l ) -gt 0 ]; then
    continue
  fi

  # Otherwise, first determine the full list of all VCFs generated for this 
  # contig in any prior joint genotyping run
  gsutil -m ls \
    $WORKSPACE_BUCKET/cromwell-execution/GnarlyJointGenotypingPart1/**/dfci-g2c.v1.$contig.*.vcf.gz \
  > $staging_dir/$contig.jg_vcfs.uris.list

  # Count the number of VCFs generated for each base workflow ID, sort s/t 
  # workflows with fewer VCFs are processed first (assuming less completion),
  # and copy all VCFs from each workflow into the staging directory
  while read wid; do
    awk -v FS="/" -v OFS="\n" -v wid="$wid" \
      '{ if ($6==wid) print $0, $0".tbi" }' \
      $staging_dir/$contig.jg_vcfs.uris.list \
    | gsutil -m cp -I \
      $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/$contig/
  done < <( awk -v FS="/" '{ print $6 }' \
            $staging_dir/$contig.jg_vcfs.uris.list \
            | sort | uniq -c | sort -nk1,1 | awk '{ print $2 }' )
done < contig_lists/dfci-g2c.v1.contigs.$WN.list

# Our next step is to find intervals covered by VCFs that successfully completed
# Gnarly joint genotyping part 1. This requires a subworkflow, below:

# Build chromosome-specific override json of VCFs and VCF indexes
while read contig; do
  # Do nothing if contig was one of the development contigs processed earlier
  if [ $( gsutil ls $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart1/$contig 2>/dev/null | wc -l ) -gt 0 ]; then
    continue
  fi

  # VCFs
  gsutil -m ls \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/$contig/**vcf.gz 2>/dev/null \
  | sort -V > $staging_dir/$contig.vcfs.list
  gsutil cp \
    $staging_dir/$contig.vcfs.list \
    $WORKSPACE_BUCKET/misc/cromwell-inputs/
done < contig_lists/dfci-g2c.v1.contigs.$WN.list

# Write template .json for input
cat << EOF > $staging_dir/UltraParallelGetVcfTerritories.inputs.template.json
{
  "UltraParallelGetVcfTerritories.g2c_analysis_docker": "vanallenlab/g2c_analysis:bd86493",
  "UltraParallelGetVcfTerritories.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "UltraParallelGetVcfTerritories.output_prefix": "dfci-g2c.v1.\$CONTIG",
  "UltraParallelGetVcfTerritories.vcf_uri_list": "$WORKSPACE_BUCKET/misc/cromwell-inputs/\$CONTIG.vcfs.list"
}
EOF

# Gather chromosomal territory covered by variant calls in finished Gnarly VCF shards
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/UltraParallelGetVcfTerritories.wdl \
  --input-json-template $staging_dir/UltraParallelGetVcfTerritories.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/GetTerritoriesGnarlyFirstPass/ \
  --contig-list <( fgrep \
                     -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
                     contig_lists/dfci-g2c.v1.contigs.$WN.list ) \
  --name GetTerritoriesGnarlyFirstPass \
  --status-tsv cromshell/progress/dfci-g2c.v1.GetTerritoriesGnarlyFirstPass.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 90 \
  --submission-gate 90 \
  --max-attempts 3

# Copy all original calling intervals to staging bucket
if ! [ -e $staging_dir/original_intervals/ ]; then
  mkdir $staging_dir/original_intervals/
fi
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/gatkhc.wgs_calling_regions.hg38.chr*.sharded.interval_list \
  $staging_dir/original_intervals/

# Build table comparing callable intervals to finished intervals for each contig
while read contig; do
  # Convert original intervals to BED
  fgrep -v "@" \
    $staging_dir/original_intervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list \
  | awk -v OFS="\t" '{ print $1, int($2), int($3) }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  > $staging_dir/original_intervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.intervals.bed

  # Get total number of callable bases
  total=$( awk '{ sum+=$3-$2 }END{ print sum }' \
             $staging_dir/original_intervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.intervals.bed )

  # Get total territory spanned by complete VCFs
  called=$( gsutil -m cat \
              $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/GetTerritoriesGnarlyFirstPass/$contig/CalcDensity/dfci-g2c.v1.$contig.density.bed.gz \
            | gunzip -c | awk '{ sum+=$3-$2 }END{ print sum }' )

  # Report
  echo -e "$contig\t$total\t$called" \
  | awk -v OFS="\t" '{ print $1, $2, $3, $3/$2 }'

done < <( fgrep -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
            contig_lists/dfci-g2c.v1.contigs.w$w.list )

# Find the list of unfinished shards per contig
# Note that this requires the prior two steps to have been run
while read contig; do
  # Only retain missed intervals ≥1kb in size
  # We assume anything smaller than this is due to boundaries at the edges of 
  # calling intervals where no variants were found
  fgrep "@" \
    $staging_dir/original_intervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list \
  > $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.preshard.interval_list
  gsutil -m cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/GetTerritoriesGnarlyFirstPass/$contig/CalcDensity/dfci-g2c.v1.$contig.density.bed.gz \
  | bedtools subtract \
    -a $staging_dir/original_intervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.intervals.bed \
    -b - \
  | awk -v OFS="\t" '{ if ($3-$2>=1000) print $1, $2, $3, "+", ". intersection ACGTmer"}' \
  >> $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.preshard.interval_list

  # If fewer than 50 shards remain, re-shard these intervals s/t >50 shards are processed.
  # This will help to solve the issues seen for LCRs clogging some shards.
  n_shards_raw=$( fgrep -v "@" $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.preshard.interval_list | wc -l )
  if [ $n_shards_raw -lt 50 ]; then
    target_bp=$( fgrep -v "@" \
                  $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.preshard.interval_list \
                | awk '{ sum+=$3-$2 }END{ print int(sum / 50) }' )
    code/scripts/split_intervals.py -v \
      -i $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.preshard.interval_list \
      -t $target_bp \
      --min-interval-size 1000 \
      -o $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.sharded.interval_list
  else
    cp \
      $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.preshard.interval_list \
      $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.sharded.interval_list
  fi
  gsutil -m cp \
    $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.sharded.interval_list \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/patch1/

done < <( fgrep -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
            contig_lists/dfci-g2c.v1.contigs.$WN.list )


####################################
# Second pass for joint genotyping #
####################################

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Reaffirm staging directory
staging_dir=staging/JGSecondPass
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write .json of contig-specific scatter counts
echo "{ " > $staging_dir/contig_variable_overrides.json
while read contig; do
  kc=$( fgrep -v "@" \
          $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.patch1.sharded.interval_list \
        | wc -l | awk '{ printf "%i\n", $1 }' )
  echo "\"$contig\" : {\"CONTIG_SCATTER_COUNT\" : $kc },"
done < <( fgrep -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
            contig_lists/dfci-g2c.v1.contigs.$WN.list ) \
| paste -s -d\  | sed 's/,$//g' \
>> $staging_dir/contig_variable_overrides.json
echo " }" >> $staging_dir/contig_variable_overrides.json

# Write template .json for input
cat << EOF > $staging_dir/GnarlyJointGenotypingPart1.patch1.inputs.template.json
{
  "GnarlyJointGenotypingPart1.callset_name": "dfci-g2c.v1.\$CONTIG.patch1",
  "GnarlyJointGenotypingPart1.dbsnp_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
  "GnarlyJointGenotypingPart1.GnarlyGenotyperFT.machine_mem_mb": 12000,
  "GnarlyJointGenotypingPart1.gnarly_scatter_count": 1,
  "GnarlyJointGenotypingPart1.import_gvcfs_batch_size": 100,
  "GnarlyJointGenotypingPart1.import_gvcfs_disk_gb": 20,
  "GnarlyJointGenotypingPart1.ImportGVCFsFT.jvm_max_mb": 25000,
  "GnarlyJointGenotypingPart1.ImportGVCFsFT.machine_mem_mb": 44000,
  "GnarlyJointGenotypingPart1.ImportGVCFsFT.n_preemptible_tries": 1,
  "GnarlyJointGenotypingPart1.intervals_already_split": true,
  "GnarlyJointGenotypingPart1.make_hard_filtered_sites": false,
  "GnarlyJointGenotypingPart1.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
  "GnarlyJointGenotypingPart1.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "GnarlyJointGenotypingPart1.ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "GnarlyJointGenotypingPart1.sample_name_map": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/dfci-g2c.v1.gatkhc.sample_map.tsv",
  "GnarlyJointGenotypingPart1.top_level_scatter_count": \$CONTIG_SCATTER_COUNT,
  "GnarlyJointGenotypingPart1.unpadded_intervals_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/patch1/gatkhc.wgs_calling_regions.hg38.\$CONTIG.patch1.sharded.interval_list"
}
EOF

# Joint genotype patch 1 per chromosome using chromsharded manager
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-hc/GnarlyJointGenotypingPart1.wdl \
  --input-json-template $staging_dir/GnarlyJointGenotypingPart1.patch1.inputs.template.json \
  --contig-variable-overrides $staging_dir/contig_variable_overrides.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotypingPatch1/ \
  --dependencies-zip gatkhc.dependencies.zip \
  --name JointGenotypingPatch1 \
  --contig-list <( fgrep -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
                     contig_lists/dfci-g2c.v1.contigs.$WN.list ) \
  --status-tsv cromshell/progress/dfci-g2c.v1.JointGenotypingPatch1.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 60 \
  --submission-gate 60 \
  --gcp-report-period 10 \
  --vm-gate 500 \
  --max-attempts 1



# ################################################
# # Clean up failed shards from joint genotyping #
# ################################################

# # In practice, it seems almost unavoidable that a handful (<5%) of shards will fail
# # either at ImportGVCFs or GnarlyGenotyper, usually due to inadequate resources.
# # Instead of increasing the resources for all tasks (and therefore increasing cost),
# # we can manually stage the VCFs per chromosome as below:

# # Note that this should be run once for each workspace for best parallelization

# # Reaffirm staging directory
# staging_dir=staging/PosthocCleanup
# if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# # Locate and stage successfully genotyped shards, then clear all prior execution buckets
# while read contig; do

#   echo $contig

#   # Find and stage successful shards
#   wid=$( tail -n1 cromshell/job_ids/dfci-g2c.v1.JointGenotyping.$contig.job_ids.list )
#   gsutil -m cp \
#     $WORKSPACE_BUCKET/cromwell-execution/GnarlyJointGenotypingPart1/$wid/**/call-GnarlyGenotyperFT/**dfci-g2c.v1.$contig.*.vcf.gz* \
#     $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/$contig/

#   # Clear all prior execution buckets
#   gsutil ls $( awk -v prefix=$WORKSPACE_BUCKET '{ print prefix"/cromwell-execution/GnarlyJointGenotypingPart1/"$1"/**" }' \
#                  cromshell/job_ids/dfci-g2c.v1.JointGenotyping.$contig.job_ids.list \
#                | paste -s - ) \
#   > uris_to_delete.list
#   cleanup_garbage

# done < <( fgrep -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
#                      contig_lists/dfci-g2c.v1.contigs.$WN.list )

# # Re-shard unsuccessful intervals and copy them to temporary bucket
# while read contig; do

#   echo $contig

#   # Reaffirm contig-specific staging directory
#   contig_dir=$staging_dir/$contig
#   if ! [ -e $contig_dir ]; then mkdir $contig_dir; fi

#   # Get expected number of shards
#   gsutil -m cp \
#     $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list \
#     $contig_dir/
#   n_shards=$( fgrep -v "@" $contig_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list | wc -l )

#   # Derive list of failed shards
#   gsutil -m ls \
#     $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/$contig/*vcf.gz \
#   | xargs -I {} basename {} \
#   | sed -e "s/\.$contig\./\t/g" -e 's/\.0\.vcf/\t/g' \
#   | awk '{ print $2 }' \
#   | sort -nk1,1 \
#   > $contig_dir/$contig.complete_shards.idx.list

#   # Infer list of shards to be rerun
#   seq 0 $((n_shards - 1)) \
#   | fgrep -xvf $contig_dir/$contig.complete_shards.idx.list \
#   > $contig_dir/$contig.failed_shards.idx.list
#   awk '{ $1+1 }' $contig_dir/$contig.failed_shards.idx.list

#   # Make list of actual intervals to be rerun
#   fgrep "@" $contig_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list \
#   > $contig_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.patch.preshard.interval_list
#   fgrep -v "@" $contig_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list \
#   | awk 'NR==FNR {a[$1+1]; next} FNR in a' \
#     $contig_dir/$contig.failed_shards.idx.list - \
#   | sort -Vk1,1 -k2,2n -k3,3n \
#   | bedtools merge -i - -o distinct -c 4,5 \
#   >> $contig_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.patch.preshard.interval_list

#   # Shard failed intervals up to a max of 5x or the contig-specific shard limit
#   n_failed_shards=$( cat $contig_dir/$contig.failed_shards.idx.list | wc -l )
#   max_permitted_shards=$n_shards
#   max_permitted_scalar=5
#   n_new_shards=$( echo -e "$(( $n_failed_shards * $max_permitted_scalar ))\n$max_permitted_shards" | sort -n | head -n1 )
#   code/scripts/split_intervals.py \
#     -i $contig_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.patch.preshard.interval_list \
#     --n-shards $n_new_shards \
#     --verbose \
#     --min-size 1000 \
#     -o $contig_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.patch.resharded.interval_list
#   gsutil cp \
#     $contig_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.patch.resharded.interval_list \
#     $WORKSPACE_BUCKET/misc/gatkhc_debug/

# done < <( fgrep -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
#                      contig_lists/dfci-g2c.v1.contigs.$WN.list )

# # Write .json of contig-specific scatter counts for second pass of joint genotyping on failed shards
# echo "{ " > $staging_dir/contig_variable_overrides.json
# while read contig; do
#   kc=$( fgrep -v "@" \
#           $staging_dir/$contig/gatkhc.wgs_calling_regions.hg38.$contig.sharded.patch.resharded.interval_list \
#         | wc -l | awk '{ printf "%i\n", $1 }' )
#   echo "\"$contig\" : {\"CONTIG_SCATTER_COUNT\" : $kc },"
# done < <( fgrep -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
#                      contig_lists/dfci-g2c.v1.contigs.$WN.list ) \
# | paste -s -d\  | sed 's/,$//g' \
# >> $staging_dir/contig_variable_overrides.json
# echo " }" >> $staging_dir/contig_variable_overrides.json

# # Write template input .json for second pass of joint genotyping on failed shards
# cat << EOF > $staging_dir/GnarlyJointGenotypingPart1Patch.inputs.template.json
# {
#   "GnarlyJointGenotypingPart1.callset_name": "dfci-g2c.v1.\$CONTIG.patch",
#   "GnarlyJointGenotypingPart1.dbsnp_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
#   "GnarlyJointGenotypingPart1.GnarlyGenotyperFT.machine_mem_mb": 16000,
#   "GnarlyJointGenotypingPart1.gnarly_scatter_count": 1,
#   "GnarlyJointGenotypingPart1.import_gvcfs_batch_size": 100,
#   "GnarlyJointGenotypingPart1.import_gvcfs_disk_gb": 20,
#   "GnarlyJointGenotypingPart1.ImportGVCFsFT.machine_mem_mb": 48000,
#   "GnarlyJointGenotypingPart1.intervals_already_split": true,
#   "GnarlyJointGenotypingPart1.make_hard_filtered_sites": false,
#   "GnarlyJointGenotypingPart1.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
#   "GnarlyJointGenotypingPart1.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
#   "GnarlyJointGenotypingPart1.ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
#   "GnarlyJointGenotypingPart1.sample_name_map": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/dfci-g2c.v1.gatkhc.sample_map.tsv",
#   "GnarlyJointGenotypingPart1.top_level_scatter_count": \$CONTIG_SCATTER_COUNT,
#   "GnarlyJointGenotypingPart1.unpadded_intervals_file": "$WORKSPACE_BUCKET/misc/gatkhc_debug/gatkhc.wgs_calling_regions.hg38.$contig.sharded.patch.resharded.interval_list"
# }
# EOF

# # Submit & monitor joint genotyping patch
# code/scripts/manage_chromshards.py \
#   --wdl code/wdl/gatk-hc/GnarlyJointGenotypingPart1.wdl \
#   --input-json-template $staging_dir/GnarlyJointGenotypingPart1Patch.inputs.template.json \
#   --contig-variable-overrides $staging_dir/contig_variable_overrides.json \
#   --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/patch/ \
#   --name JointGenotypingPatch \
#   --contig-list <( fgrep -xvf contig_lists/dfci-g2c.v1.contigs.dev.list \
#                               contig_lists/dfci-g2c.v1.contigs.$WN.list ) \
#   --status-tsv cromshell/progress/dfci-g2c.v1.JointGenotypingPatch.progress.tsv \
#   --workflow-id-log-prefix "dfci-g2c.v1" \
#   --outer-gate 60 \
#   --submission-gate 30 \
#   --vm-gate 1000 \
#   --max-attempts 3

# # Relocate the outputs from successful patch jobs to the main staging bucket
# # TODO: implement this

# # TODO: probably will need to handle shards that continue to fail
# # TBD on how much to automate this

###############
# VCF cleanup #
###############

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Refresh staging directory
staging_dir=staging/PosthocCleanup
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Build chromosome-specific override json of VCFs and VCF indexes
echo "{}" > $staging_dir/contig_variable_overrides.json
while read contig; do
  # VCFs
  gsutil -m ls \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/$contig/**vcf.gz \
  | sort -V > $staging_dir/$contig.vcfs.list

  # VCF indexes
  awk '{ print $1".tbi" }' $staging_dir/$contig.vcfs.list \
  > $staging_dir/$contig.vcf_idxs.list

  # Write .json snippet for variable overrides for this contig
  cat << EOF > $staging_dir/$contig.overrides.json
{
  "$contig" : {
      "CONTIG_VCFS": $( collapse_txt $staging_dir/$contig.vcfs.list ),
      "CONTIG_VCF_IDXS": $( collapse_txt $staging_dir/$contig.vcf_idxs.list )
    }
}
EOF
  
  # Update main .json
  code/scripts/update_json.py \
    -i $staging_dir/contig_variable_overrides.json \
    -u $staging_dir/$contig.overrides.json \
    -o $staging_dir/contig_variable_overrides.json
done < contig_lists/dfci-g2c.v1.contigs.$WN.list

# Write template .json for input
cat << EOF > $staging_dir/PosthocCleanupPart1.inputs.template.json
{
  "PosthocCleanupPart1.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PosthocCleanupPart1.g2c_pipeline_docker": "vanallenlab/g2c_pipeline:sv_counting",
  "PosthocCleanupPart1.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "PosthocCleanupPart1.NormalizeVcf.mem_gb": 31,
  "PosthocCleanupPart1.output_prefix": "dfci-g2c.v1.\$CONTIG",
  "PosthocCleanupPart1.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "PosthocCleanupPart1.vcfs": \$CONTIG_VCFS,
  "PosthocCleanupPart1.vcf_idxs": \$CONTIG_VCF_IDXS
}
EOF

# Perform post hoc VCF cleanup (split multiallelics, minimize indel representation)
# Also count variants by type per sample (needed for outlier definition; see below)
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-hc/PosthocCleanupPart1.wdl \
  --input-json-template $staging_dir/PosthocCleanupPart1.inputs.template.json \
  --contig-variable-overrides $staging_dir/contig_variable_overrides.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart1/ \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.PosthocCleanupPart1.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 45 \
  --max-attempts 2


###########################
# Outlier sample analysis #
###########################

# Note that this section only needs to be run from one workspace for the entire cohort

# Reaffirm staging directory
staging_dir=staging/PosthocCleanup
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Sum variant counts for all contigs
for k in $( seq 1 22 ) X Y; do
  contig="chr$k"
  gsutil cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart1/$contig/PosthocCleanupPart1.$contig.outputs.json \
  | jq .counts_per_sample | tr -d '"'
done | gsutil -m cp -I $staging_dir/
code/scripts/sum_svcounts.py \
  --outfile $staging_dir/dfci-g2c.v1.gatkhc.PosthocCleanupPart1.counts.tsv \
  $staging_dir/dfci-g2c.v1.chr*.norm.counts.tsv

# Ensure R packages are installed
. code/refs/install_packages.sh R

# Copy metadata from end of GATK-SV pipeline
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
  ./

# Regenerate ancestry labels split by cohort
pop_idx=$( zcat dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
           | sed -n '1p' | sed 's/\t/\n/g' \
           | awk '{ if ($1=="intake_qc_pop") print NR }' )
zcat dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz | sed '1d' \
| awk -v idx=$pop_idx -v FS="\t" -v OFS="\t" '{ if ($3!="aou") $3="oth"; print $1, $idx"_"$3 }' \
| cat <( echo -e "sample_id\tlabel" ) - \
> $staging_dir/dfci-g2c.intake_pop_labels.aou_split.tsv

# Define outliers
code/scripts/define_variant_count_outlier_samples.R \
  --counts-tsv $staging_dir/dfci-g2c.v1.gatkhc.PosthocCleanupPart1.counts.tsv \
  --sample-labels-tsv $staging_dir/dfci-g2c.intake_pop_labels.aou_split.tsv \
  --n-iqr 4 \
  --plot \
  --plot-title-prefix "GATK" \
  --out-prefix $staging_dir/dfci-g2c.v1.gatkhc.posthoc_outliers

# Get list of samples that were considered for joint genotyping
gsutil cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/dfci-g2c.v1.gatkhc.sample_map.tsv \
| cut -f1 | sort -V | uniq \
> $staging_dir/dfci-g2c.v1.gatkhc.samples.list

# Add non-technical samples to exclude due to age data becoming available
age_cidx=$( zcat dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
            | head -n1 | sed 's/\t/\n/g' | awk '{ if ($1=="age") print NR }'  )
zcat dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
| awk -v FS="\t" -v cidx=$age_cidx '{ if ($cidx < 18) print $1 }' \
| fgrep -wf $staging_dir/dfci-g2c.v1.gatkhc.samples.list \
| cat - $staging_dir/dfci-g2c.v1.gatkhc.posthoc_outliers.outliers.samples.list \
| sort -V | uniq \
> $staging_dir/dfci-g2c.v1.gatkhc.posthoc_outliers.outliers.samples.list2
mv $staging_dir/dfci-g2c.v1.gatkhc.posthoc_outliers.outliers.samples.list2 \
  $staging_dir/dfci-g2c.v1.gatkhc.posthoc_outliers.outliers.samples.list

# Update sample metadata with posthoc outlier failure labels
code/scripts/append_qc_fail_metadata.R \
  --qc-tsv dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
  --new-column-name gatkhc_posthoc_qc_pass \
  --all-samples-list $staging_dir/dfci-g2c.v1.gatkhc.samples.list \
  --fail-samples-list $staging_dir/dfci-g2c.v1.gatkhc.posthoc_outliers.outliers.samples.list \
  --outfile dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv
gzip -f dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv

# Compress and archive outlier data for future reference
cd $staging_dir && \
tar -czvf dfci-g2c.v1.gatkhc.posthoc_outliers.tar.gz dfci-g2c.v1.gatkhc.posthoc_outliers* && \
gsutil -m cp \
  dfci-g2c.v1.gatkhc.posthoc_outliers.tar.gz \
  dfci-g2c.v1.gatkhc.posthoc_outliers.outliers.samples.list \
  ~/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/ && \
cd ~

# Replot sample QC after excluding outliers above
qcplotdir=dfci-g2c.phase1.gatkhc_posthoc_qc_pass.plots
if [ ! -e $qcplotdir ]; then mkdir $qcplotdir; fi
code/scripts/plot_intake_qc.R \
  --qc-tsv dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --pass-column filtersites_qc_pass \
  --pass-column gatksv_posthoc_qc_pass \
  --pass-column gatkhc_posthoc_qc_pass \
  --out-prefix $qcplotdir/dfci-g2c.phase1.gatkhc_posthoc_qc_pass
tar -czvf dfci-g2c.phase1.gatkhc_posthoc_qc_pass.plots.tar.gz $qcplotdir
gsutil -m cp \
  dfci-g2c.phase1.gatkhc_posthoc_qc_pass.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/gatkhc_qc/



#############################################################
# Exclude outlier samples and apply site-level hard filters #
#############################################################

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Reaffirm staging directory
staging_dir=staging/PosthocCleanup
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Build chromosome-specific override json of VCFs and VCF indexes
add_contig_vcfs_to_chromshard_overrides_json \
  $staging_dir/PosthocCleanupPart2.contig_variable_overrides.json \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart1 \
  normalized_vcfs \
  normalized_vcf_idxs

# Write template input .json for outlier exclusion task
cat << EOF > $staging_dir/PosthocCleanupPart2.inputs.template.json
{
  "PosthocCleanupPart2.CleanupPart2.mem_gb": 7.5,
  "PosthocCleanupPart2.CleanupPart2.n_cpu": 4,
  "PosthocCleanupPart2.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PosthocCleanupPart2.exclude_samples_list": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.v1.gatkhc.posthoc_outliers.outliers.samples.list",
  "PosthocCleanupPart2.vcfs": \$CONTIG_VCFS,
  "PosthocCleanupPart2.vcf_idxs": \$CONTIG_VCF_IDXS
}
EOF

# Submit outlier exclusion & hard filter task using chromsharded manager
# Reminder that this manager script handles staging & cleanup too
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-hc/PosthocCleanupPart2.wdl \
  --input-json-template $staging_dir/PosthocCleanupPart2.inputs.template.json \
  --contig-variable-overrides $staging_dir/PosthocCleanupPart2.contig_variable_overrides.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart2/ \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.PosthocCleanupPart2.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --max-attempts 2


###########################################
# Exclude outliers also from GATK-SV VCFs #
###########################################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Reaffirm staging directory
staging_dir=staging/PosthocCleanup
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for outlier exclusion & hard filter task
cat << EOF > $staging_dir/ExcludeSnvOutliersFromSvCallset.inputs.template.json
{
  "PosthocHardFilterPart2.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PosthocHardFilterPart2.exclude_samples_list": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.v1.gatkhc.posthoc_outliers.outliers.samples.list",
  "PosthocHardFilterPart2.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/CollapseRedundantSvs/\$CONTIG/RC3/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.identical.reclustered.vcf.gz",
  "PosthocHardFilterPart2.vcf_idx": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/CollapseRedundantSvs/\$CONTIG/RC3/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.identical.reclustered.vcf.gz.tbi"
}
EOF

# Submit outlier exclusion & hard filter task using chromsharded manager
# Reminder that this manager script handles staging & cleanup too
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/PosthocHardFilterPart2.wdl \
  --input-json-template $staging_dir/ExcludeSnvOutliersFromSvCallset.inputs.template.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset \
  --name ExcludeSnvOutliersFromSvCallset \
  --status-tsv cromshell/progress/dfci-g2c.v1.ExcludeSnvOutliersFromSvCallset.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --max-attempts 2 


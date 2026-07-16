#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-SV cohort mode pipeline on G2C phase 1

# Note that this code is designed to be run inside the AoU Researcher Workbench
# See gatksv_bash_utils.sh for custom function definitions used below


# Note 2: this code was defined for the AoU RW v1.0, and it may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in cromshell cromshell/inputs cromshell/job_ids cromshell/progress; do
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
. code/refs/gatksv_bash_utils.sh

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Create dependencies .zip for all GATK-SV module submissions
cd code/wdl/gatk-sv && \
zip gatksv.dependencies.zip *.wdl && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for G2C-specific workflow submissions
cd code/wdl/pancan_germline_wgs && \
zip g2c.dependencies.zip *.wdl && \
mv g2c.dependencies.zip ~/ && \
cd ~

# Ensure Cromwell/Cromshell are configured
code/scripts/setup_cromshell.py

# Install necessary packages
. code/refs/install_packages.sh python

# Copy sample & batch information
gsutil -m cp -r \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/batch_info \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
  ~/


##################
# 03 | TrainGCNV #
##################

module_submission_routine_all_batches 03


############################
# 04 | GatherBatchEvidence #
############################

module_submission_routine_all_batches 04


#####################
# 05 | ClusterBatch #
#####################

module_submission_routine_all_batches 05


##################################
# 05B | ExcludeClusteredOutliers #
##################################

# Note: this is not a canonical GATK-SV module and was instituted specifically for G2C

# First, we need to define a list of outlier samples to be excluded
# This needs to only be run once from one of the AoU workspaces

# Install R packages
. code/refs/install_packages.sh R

# Collect SV count data for each sample per algorithm
if ! [ -e data/gatksv_05B_outliers ]; then
  mkdir gatksv_05B_outliers
fi
for alg in depth manta melt wham; do
  gsutil -m cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/05/*/PlotSVCountsPerSample/CountSVsPerSamplePerType/**.cluster_batch.$alg.svcounts.txt \
  | sort -Vrk1,1 -k2,2V | uniq \
  > gatksv_05B_outliers/dfci-g2c.05_sv_counts.$alg.tsv
done

# Define outlier samples cohort-wide as: 
# * Q3 + 6 IQR for manta, wham, or melt
# * Q3 + 10 IQR for depth CNVs
pop_idx=$( zcat dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
           | sed -n '1p' | sed 's/\t/\n/g' \
           | awk '{ if ($1=="intake_qc_pop") print NR }' )
zcat dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz | sed '1d' \
| awk -v idx=$pop_idx -v FS="\t" -v OFS="\t" '{ if ($3!="aou") $3="oth"; print $1, $idx"_"$3 }' \
| cat <( echo -e "sample_id\tlabel" ) - \
> gatksv_05B_outliers/dfci-g2c.intake_pop_labels.aou_split.tsv
for alg in depth manta melt wham; do
  niqr=6
  case $alg in
    "depth")
      niqr=10
      prefix="Depth"
      ;;
    "manta")
      prefix="Manta"
      ;;
    "melt")
      prefix="Melt"
      ;;
    "wham")
      prefix="Wham"
      ;;
  esac
  code/scripts/define_variant_count_outlier_samples.R \
    --counts-tsv gatksv_05B_outliers/dfci-g2c.05_sv_counts.$alg.tsv \
    --sample-labels-tsv gatksv_05B_outliers/dfci-g2c.intake_pop_labels.aou_split.tsv \
    --n-iqr $niqr \
    --no-lower-filter \
    --plot \
    --plot-title-prefix $prefix \
    --out-prefix gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.$alg
done

# Merge outlier sample lists
cat gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.*.outliers.samples.list \
| sort -V | uniq \
> gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.all_outlier_samples.list

# Update sample metadata with 05B outlier failure labels
code/scripts/append_qc_fail_metadata.R \
  --qc-tsv dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
  --new-column-name clusterbatch_qc_pass \
  --all-samples-list <( cat batch_info/sample_lists/* ) \
  --fail-samples-list gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.all_outlier_samples.list \
  --outfile dfci-g2c.sample_meta.post_clusterbatch.tsv
gzip -f dfci-g2c.sample_meta.post_clusterbatch.tsv

# Compress and archive outlier data for future reference
tar -czvf gatksv_05B_outliers.tar.gz gatksv_05B_outliers
gsutil -m cp \
  gatksv_05B_outliers.tar.gz \
  gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.all_outlier_samples.list \
  dfci-g2c.sample_meta.post_clusterbatch.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/

# Replot sample QC after excluding outliers above
qcplotdir=dfci-g2c.phase1.clusterbatch_qc_pass.plots
if [ ! -e $qcplotdir ]; then mkdir $qcplotdir; fi
code/scripts/plot_intake_qc.R \
  --qc-tsv dfci-g2c.sample_meta.post_clusterbatch.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --out-prefix $qcplotdir/dfci-g2c.phase1.clusterbatch_qc_pass
tar -czvf dfci-g2c.phase1.clusterbatch_qc_pass.plots.tar.gz $qcplotdir
gsutil -m cp \
  dfci-g2c.phase1.clusterbatch_qc_pass.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/gatksv_qc/

# Launch workflows to remove outlier samples from clustered VCFs for each batch
module_submission_routine_all_batches 05B


########################
# 05C | ReclusterBatch #
########################

# Note: this is not a canonical GATK-SV module and was instituted specifically for G2C

module_submission_routine_all_batches 05C


#############################
# 06 | GenerateBatchMetrics #
#############################

module_submission_routine_all_batches 06


#########################
# 07 | FilterBatchSites #
#########################

module_submission_routine_all_batches 07


###########################
# 08 | FilterBatchSamples #
###########################

# Note: similar to 05B above, this step uses custom outlier definitions
# The below code needs to be run just once in a single workspace

# This also requires all packages installed and outputs as for 05B above
# If necessary, rerun the 05B library installation steps before proceeding

# Collect SV count data for each sample per algorithm
if ! [ -e data/gatksv_08_outliers ]; then
  mkdir gatksv_08_outliers
fi
for alg in depth manta melt wham; do
  gsutil -m cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/07/*/PlotSVCountsPerSample/CountSVsPerSamplePerType/**.$alg.with_evidence.svcounts.txt \
  | sort -Vrk1,1 -k2,2V | uniq \
  > gatksv_08_outliers/dfci-g2c.08_sv_counts.$alg.tsv
done

# Ensure R packages are installed
. code/refs/install_packages.sh R

# Define outlier samples cohort-wide as Q3 + 6 IQR for all algorithms
for alg in depth manta melt wham; do
  case $alg in
    "depth")
      prefix="Depth"
      ;;
    "manta")
      prefix="Manta"
      ;;
    "melt")
      prefix="Melt"
      ;;
    "wham")
      prefix="Wham"
      ;;
  esac
  code/scripts/define_variant_count_outlier_samples.R \
    --counts-tsv gatksv_08_outliers/dfci-g2c.08_sv_counts.$alg.tsv \
    --sample-labels-tsv gatksv_05B_outliers/dfci-g2c.intake_pop_labels.aou_split.tsv \
    --n-iqr 6 \
    --no-lower-filter \
    --plot \
    --plot-title-prefix $prefix \
    --out-prefix gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.$alg
done

# Merge outlier sample lists
cat gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.*.outliers.samples.list \
| sort -V | uniq \
> gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.all_outlier_samples.list

# Update sample metadata with 08 outlier failure labels
code/scripts/append_qc_fail_metadata.R \
  --qc-tsv dfci-g2c.sample_meta.post_clusterbatch.tsv.gz \
  --new-column-name filtersites_qc_pass \
  --all-samples-list <( cat batch_info/sample_lists/* | fgrep -wvf \
                        gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.all_outlier_samples.list ) \
  --fail-samples-list gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.all_outlier_samples.list \
  --outfile dfci-g2c.sample_meta.post_filtersites.tsv
gzip -f dfci-g2c.sample_meta.post_filtersites.tsv

# Compress and archive outlier data for future reference
tar -czvf gatksv_08_outliers.tar.gz gatksv_08_outliers
gsutil -m cp \
  gatksv_08_outliers.tar.gz \
  gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.all_outlier_samples.list \
  dfci-g2c.sample_meta.post_filtersites.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/

# Replot sample QC after excluding outliers above
qcplotdir=dfci-g2c.phase1.filtersites_qc_pass.plots
if [ ! -e $qcplotdir ]; then mkdir $qcplotdir; fi
code/scripts/plot_intake_qc.R \
  --qc-tsv dfci-g2c.sample_meta.post_filtersites.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --pass-column filtersites_qc_pass \
  --out-prefix $qcplotdir/dfci-g2c.phase1.filtersites_qc_pass
tar -czvf dfci-g2c.phase1.filtersites_qc_pass.plots.tar.gz $qcplotdir
gsutil -m cp \
  dfci-g2c.phase1.filtersites_qc_pass.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/gatksv_qc/

# Launch workflows to remove outlier samples from filtered VCFs for each batch
module_submission_routine_all_batches 08


########################
# 09 | MergeBatchSites #
########################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 09

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.09-MergeBatchSites.job_ids.list )

# Once complete, stage outputs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.09-MergeBatchSites.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/09/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell*/MergeBatchSites/** >> uris_to_delete.list
cleanup_garbage


######################
# 10 | GenotypeBatch #
######################

module_submission_routine_all_batches 10


#######################
# 11 | RegenotypeCNVs #
#######################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 11

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.11-RegenotypeCNVs.job_ids.list )

# Once complete, stage outputs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.11-RegenotypeCNVs.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/11/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell*/RegenotypeCNVs/** >> uris_to_delete.list
cleanup_garbage


#######################
# 12 | CombineBatches #
#######################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Note 2: this module is handled differently by submit_cohort_module since it's
# parallelized by chromosome with 24 independent submissions

# All cleanup and tracking is handled by a helper routine within submit_cohort_module

submit_cohort_module 12


###############################
# 13 | ResolveComplexVariants #
###############################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 13

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.13-ResolveComplexVariants.job_ids.list )

# Once complete, stage outputs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.13-ResolveComplexVariants.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/13/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell*/ResolveComplexVariants/** >> uris_to_delete.list
cleanup_garbage


###############################
# 14A | FilterCoverageSamples #
###############################

# Note: this is not a canonical GATK-SV module and was instituted specifically for G2C
# All this does is remove outliers from the median coverage file to avoid 
# RDtest normalization errors in module 14

# Get list of samples present in VCF at this stage
if [ ! -e staging/14A-FilterCoverageSamples ]; then
  mkdir staging/14A-FilterCoverageSamples
fi
gsutil cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/13/dfci-g2c.v1.reshard_vcf.chrY.resharded.vcf.gz \
| gunzip -c | head -n300 | fgrep -v "##" | fgrep "#" | cut -f10- \
| sed 's/\t/\n/g' | sort -V | uniq \
> staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list
gsutil -m cp \
  staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/

module_submission_routine_all_batches 14A


################################
# 14 | GenotypeComplexVariants #
################################

# Note: this module only needs to be run once in one workspace for the whole cohort

# # All cleanup and tracking is handled by a helper routine within submit_cohort_module

# Subset .ped file to those present in VCF
samples14A_list=staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list
if [ ! -e $samples14A_list ]; then
  echo "ERROR: cannot locate $samples14A_list. Regenerate if necessary."
else
  if [ ! -e staging/14-GenotypeComplexVariants ]; then
    mkdir staging/14-GenotypeComplexVariants
  fi
  gsutil -m cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped \
  | fgrep -wf staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list \
  > staging/14-GenotypeComplexVariants/dfci-g2c.all_samples.gatksv_module14.ped
  gsutil cp \
    staging/14-GenotypeComplexVariants/dfci-g2c.all_samples.gatksv_module14.ped \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/
fi

# Submit workflow
submit_cohort_module 14

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.14-GenotypeComplexVariants.job_ids.list )

# Once complete, stage outputs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.14-GenotypeComplexVariants.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/14/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell*/GenotypeComplexVariants/** >> uris_to_delete.list
cleanup_garbage


#################
# 15 | CleanVcf #
#################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 15

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.15-CleanVcf.job_ids.list )

# CleanVcf processes all chromosomes in parallel but is designed to concatenate 
# all chromosomes into a single output VCF. However, since all of our downstream
# steps are parallel by chromosome, there is no point in waiting for CleanVcf to
# finish concatenating the massive overall VCF. Thus, we can determine whether 
# CleanVcf is finished by checking if the concatenation task is present in the
# execution bucket, in which case we can kill the job and relocalize all of the
# chromosome-sharded outputs
cvcf_exec_base=$WORKSPACE_BUCKET/cromwell-execution/CleanVcf/$( tail -n1 cromshell/job_ids/dfci-g2c.v1.15-CleanVcf.job_ids.list )
if [ $( gsutil ls $cvcf_exec_base/call-ConcatCleanedVcfs | wc -l ) -gt 0 ]; then
  for k in $( seq 0 23 ); do
    contig_wid=$( basename $( gsutil ls $cvcf_exec_base/call-CleanVcfChromosome/shard-$k/CleanVcfChromosome/ ) )
    cromshell -t 120 --no_turtle -mc list-outputs $contig_wid \
    | awk '{ print $2 }' | gsutil -m cp -I \
      $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/15/
  done
  cromshell -t 120 abort $( tail -n1 cromshell/job_ids/dfci-g2c.v1.15-CleanVcf.job_ids.list )
fi

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell*/CleanVcf/** >> uris_to_delete.list
cleanup_garbage


##############################
# 16 | RefineComplexVariants #
##############################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Subset batch sample membership files to those present in VCF
samples14A_list=staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list
if [ ! -e $samples14A_list ]; then
  echo "ERROR: cannot locate $samples14A_list. Regenerate if necessary."
else
  for dir in 16-RefineComplexVariants \
             16-RefineComplexVariants/module14A_batch_sample_lists; do
    if [ ! -e staging/$dir ]; then mkdir staging/$dir; fi
  done
  while read bid; do
    fgrep -xf $samples14A_list batch_info/sample_lists/$bid.samples.list \
    > staging/16-RefineComplexVariants/module14A_batch_sample_lists/$bid.gatksv_module14A.samples.list
  done < batch_info/dfci-g2c.gatk-sv.batches.list
  gsutil cp -r \
    staging/16-RefineComplexVariants/module14A_batch_sample_lists \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/
fi

# Note 2: this module is handled differently by submit_cohort_module since it's
# parallelized by chromosome with 24 independent submissions

# All cleanup and tracking is handled by a helper routine within submit_cohort_module

submit_cohort_module 16


###########################################################
# Post hoc site hard filters and outlier sample exclusion #
###########################################################

# Note that this section only needs to be run from one workspace for the entire cohort

# Write template input .json for hard filters, part 1
staging_dir=staging/posthoc_filter
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi
mkdir $staging_dir
cat << EOF > $staging_dir/PosthocHardFilterPart1.inputs.template.json
{
  "PosthocHardFilterPart1.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PosthocHardFilterPart1.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/16/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.cpx_refined.vcf.gz",
  "PosthocHardFilterPart1.vcf_idx": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/16/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.cpx_refined.vcf.gz.tbi"
}
EOF

# Submit first round of site hard filters using chromsharded manager
# Reminder that this manager script handles staging & cleanup too
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/PosthocHardFilterPart1.wdl \
  --input-json-template $staging_dir/PosthocHardFilterPart1.inputs.template.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/PosthocHardFilterPart1 \
  --name PosthocHardFilterPart1 \
  --status-tsv cromshell/progress/dfci-g2c.v1.PosthocHardFilterPart1.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --max-attempts 3

# Write input .json for SV counting task
for k in $( seq 1 22 ) X Y; do
  gsutil cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/PosthocHardFilterPart1/chr$k/PosthocHardFilterPart1.chr$k.outputs.json \
  | jq .filtered_vcf | tr -d '"' \
  >> $staging_dir/vcfs.list
  gsutil cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/PosthocHardFilterPart1/chr$k/PosthocHardFilterPart1.chr$k.outputs.json \
  | jq .filtered_vcf_idx | tr -d '"' \
  >> $staging_dir/vcf_idxs.list
done
cat << EOF > cromshell/inputs/count_svs_posthoc.inputs.json
{
  "CountSvsPerSample.ShardVcf.disk_gb": 100,
  "CountSvsPerSample.ShardVcf.n_preemptible": 0,
  "CountSvsPerSample.g2c_pipeline_docker": "vanallenlab/g2c_pipeline:sv_counting",
  "CountSvsPerSample.output_prefix": "dfci-g2c.v1.gatksv_postCleanVcf",
  "CountSvsPerSample.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-01-14-v1.0.1-88dbd052",
  "CountSvsPerSample.vcfs": $( collapse_txt $staging_dir/vcfs.list ),
  "CountSvsPerSample.vcf_idxs": $( collapse_txt $staging_dir/vcf_idxs.list )
}
EOF

# Submit SV counting task
cromshell --no_turtle -t 120 -mc submit \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip g2c.dependencies.zip \
  code/wdl/gatk-sv/CountSvsPerSample.wdl \
  cromshell/inputs/count_svs_posthoc.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/count_svs_posthoc.job_ids.list

# Monitor SV counting task
monitor_workflow $( tail -n1 cromshell/job_ids/count_svs_posthoc.job_ids.list )

# Ensure R packages are installed
. code/refs/install_packages.sh R

# Once complete, download SV counts per sample, collapse CPX+INV, and exclude CTX and mCNVs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/count_svs_posthoc.job_ids.list ) \
| awk '{ print $NF }' | gsutil cp -I $staging_dir/
sed 's/\tINV\t\|\tCPX\t/\tINV_CPX\t/g' \
  $staging_dir/dfci-g2c.v1.gatksv_postCleanVcf.counts.tsv \
| awk -v FS="\t" -v OFS="\t" '{ if ($2 !~ /CTX|BND|CNV/) print }' \
> $staging_dir/dfci-g2c.v1.gatksv_postCleanVcf.counts.subsetted.tsv
code/scripts/sum_svcounts.py \
  --outfile $staging_dir/dfci-g2c.v1.gatksv_postCleanVcf.counts.subsetted.collapsed.tsv \
  $staging_dir/dfci-g2c.v1.gatksv_postCleanVcf.counts.subsetted.tsv

# Regenerate ancestry labels split by cohort
pop_idx=$( zcat dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
           | sed -n '1p' | sed 's/\t/\n/g' \
           | awk '{ if ($1=="intake_qc_pop") print NR }' )
zcat dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz | sed '1d' \
| awk -v idx=$pop_idx -v FS="\t" -v OFS="\t" '{ if ($3!="aou") $3="oth"; print $1, $idx"_"$3 }' \
| cat <( echo -e "sample_id\tlabel" ) - \
> $staging_dir/dfci-g2c.intake_pop_labels.aou_split.tsv

# Define outliers
code/scripts/define_variant_count_outlier_samples.R \
  --counts-tsv $staging_dir/dfci-g2c.v1.gatksv_postCleanVcf.counts.subsetted.collapsed.tsv \
  --sample-labels-tsv $staging_dir/dfci-g2c.intake_pop_labels.aou_split.tsv \
  --n-iqr 4 \
  --plot \
  --plot-title-prefix "GATKSV" \
  --out-prefix $staging_dir/dfci-g2c.v1.gatksv.posthoc_outliers

# Update sample metadata with posthoc outlier failure labels
code/scripts/append_qc_fail_metadata.R \
  --qc-tsv dfci-g2c.sample_meta.post_filtersites.tsv.gz \
  --new-column-name gatksv_posthoc_qc_pass \
  --all-samples-list <( cat batch_info/sample_lists/* | fgrep -wf \
                        staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list ) \
  --fail-samples-list $staging_dir/dfci-g2c.v1.gatksv.posthoc_outliers.outliers.samples.list \
  --outfile dfci-g2c.sample_meta.posthoc_outliers.tsv

# Update CEPH phenotypes in sample metadata .tsv
# This is necessary because we gained access to the CEPH cancer data now,
# chronologically months after sample processing had started. There will be a slight
# discrepancy between batching case:control labels for the 130 CEPH cancer cases
gsutil cp \
  gs://dfci-g2c-inputs/phenotypes/ceph.phenos.tsv.gz \
  ./
code/scripts/update_metadata_phenotypes.R \
  --metadata-tsv dfci-g2c.sample_meta.posthoc_outliers.tsv \
  --phenotypes-tsv ceph.phenos.tsv.gz \
  --out-tsv dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv
gzip -f dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv

# Compress and archive outlier data for future reference
cd $staging_dir && \
tar -czvf dfci-g2c.v1.gatksv.posthoc_outliers.tar.gz dfci-g2c.v1.gatksv.posthoc_outliers* && \
gsutil -m cp \
  dfci-g2c.v1.gatksv.posthoc_outliers.tar.gz \
  dfci-g2c.v1.gatksv.posthoc_outliers.outliers.samples.list \
  ~/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/ && \
cd ~

# Replot sample QC after excluding outliers above
qcplotdir=dfci-g2c.phase1.gatksv_posthoc_qc_pass.plots
if [ ! -e $qcplotdir ]; then mkdir $qcplotdir; fi
code/scripts/plot_intake_qc.R \
  --qc-tsv dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --pass-column filtersites_qc_pass \
  --pass-column gatksv_posthoc_qc_pass \
  --out-prefix $qcplotdir/dfci-g2c.phase1.gatksv_posthoc_qc_pass
tar -czvf dfci-g2c.phase1.gatksv_posthoc_qc_pass.plots.tar.gz $qcplotdir
gsutil -m cp \
  dfci-g2c.phase1.gatksv_posthoc_qc_pass.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/gatksv_qc/

# Write template input .json for outlier exclusion & hard filter task
cat << EOF > $staging_dir/PosthocHardFilterPart2.inputs.template.json
{
  "PosthocHardFilterPart2.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PosthocHardFilterPart2.exclude_samples_list": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/dfci-g2c.v1.gatksv.posthoc_outliers.outliers.samples.list",
  "PosthocHardFilterPart2.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/PosthocHardFilterPart1/\$CONTIG/HardFilterPart1/dfci-g2c.v1.\$CONTIG.cpx_refined.posthoc_filtered.vcf.gz",
  "PosthocHardFilterPart2.vcf_idx": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/PosthocHardFilterPart1/\$CONTIG/HardFilterPart1/dfci-g2c.v1.\$CONTIG.cpx_refined.posthoc_filtered.vcf.gz.tbi"
}
EOF

# Submit outlier exclusion & hard filter task using chromsharded manager
# Reminder that this manager script handles staging & cleanup too
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/PosthocHardFilterPart2.wdl \
  --input-json-template $staging_dir/PosthocHardFilterPart2.inputs.template.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/PosthocHardFilterPart2 \
  --name PosthocHardFilterPart2 \
  --status-tsv cromshell/progress/dfci-g2c.v1.PosthocHardFilterPart2.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --max-attempts 3


#####################
# 17 | JoinRawCalls #
#####################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 17

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.17-JoinRawCalls.job_ids.list ) \
  20

# Once complete, stage outputs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.17-JoinRawCalls.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/17/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell*/JoinRawCalls/** >> uris_to_delete.list
cleanup_garbage


######################
# 18 | SVConcordance #
######################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Note 2: this module is handled differently by submit_cohort_module since it's
# parallelized by chromosome with 24 independent submissions

# All cleanup and tracking is handled by a helper routine within submit_cohort_module

submit_cohort_module 18


########################
# 19 | FilterGenotypes #
########################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Note 2: this module is handled differently by submit_cohort_module since it's
# parallelized by chromosome with 24 independent submissions

# All cleanup and tracking is handled by a helper routine within submit_cohort_module

submit_cohort_module 19



#####################################
# Collapse quasi-redundant variants #
#####################################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Write template input .json for reclustering
staging_dir=staging/posthoc_recluster
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi
cat << EOF > $staging_dir/CollapseRedundantSvs.inputs.template.json
{
  "CollapseRedundantSvs.g2c_analysis_docker": "vanallenlab/g2c_analysis:e246265",
  "CollapseRedundantSvs.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/19/\$CONTIG/RecalibrateGq/ConcatVcfs/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.vcf.gz",
  "CollapseRedundantSvs.vcf_idx": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/19/\$CONTIG/RecalibrateGq/ConcatVcfs/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.vcf.gz.tbi"
}
EOF

# Submit, monitor, and stage/cleanup redundant variant reclustering
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/CollapseRedundantSvs.wdl \
  --input-json-template $staging_dir/CollapseRedundantSvs.inputs.template.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/CollapseRedundantSvs \
  --name CollapseRedundantSvs \
  --status-tsv cromshell/progress/dfci-g2c.v1.CollapseRedundantSvs.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --submission-gate 0 \
  --max-attempts 3

# Confirm that the reclustering procedure didn't result in major VCF changes
for k in $( seq 1 22 ) X Y; do
  for wrapper in 1; do
    CONTIG="chr$k"
    echo $CONTIG
    gsutil -m cat \
      $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/19/$CONTIG/RecalibrateGq/ConcatVcfs/dfci-g2c.v1.$CONTIG.concordance.gq_recalibrated.vcf.gz \
    | bcftools query -f '%ID\n' | wc -l
    gsutil -m cat \
      $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/CollapseRedundantSvs/$CONTIG/RC3/dfci-g2c.v1.$CONTIG.concordance.gq_recalibrated.identical.reclustered.vcf.gz \
    | bcftools query -f '%ID\n' | wc -l
  done | paste - - - | awk -v OFS="\t" '{ print $1, $2, $3, $3-$2, 100*(($3-$2)/$2)"%" }'
done

# Test downstream tool compatability by trying to annotate all CPX variants

# Make tarball with SV annotation dependencies
cd code/wdl/gatk-sv/svannotation_v1.1 && \
zip svannotation.dependencies.zip *.wdl && \
mv svannotation.dependencies.zip ~/ && \
cd ~

# Extract CPX variants from chr1 (we know there's a INVdup at the end that has been problematic)
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/CollapseRedundantSvs/chr1/RC3/dfci-g2c.v1.chr1.concordance.gq_recalibrated.identical.reclustered.vcf.gz \
| bcftools view \
  -i 'INFO/SVTYPE = "CPX"' \
  -Oz -o scratch/dfci-g2c.v1.sv.redundant_collapsed.cpx_only.chr1.vcf.gz
tabix -p vcf -f scratch/dfci-g2c.v1.sv.redundant_collapsed.cpx_only.chr1.vcf.gz
gsutil -m cp \
  scratch/dfci-g2c.v1.sv.redundant_collapsed.cpx_only.chr1.vcf.gz* \
  $WORKSPACE_BUCKET/scratch/

# Write .json input for chr1 annotation test
cat << EOF > cromshell/inputs/AnnotateVcf.inputs.chr1_test.json
{
  "AnnotateVcf.contig_list": "gs://dfci-g2c-refs/hg38/contig_lists/chr1.list",
  "AnnotateVcf.external_af_population": ["ALL","AFR","AMR","EAS","EUR","MID","FIN","ASJ","RMI","SAS","AMI"],
  "AnnotateVcf.external_af_ref_bed": "gs://gatk-sv-resources-public/gnomad_AF/gnomad_v4_SV.Freq.tsv.gz",
  "AnnotateVcf.external_af_ref_prefix": "gnomad_v4.1_sv",
  "AnnotateVcf.gatk_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT",
  "AnnotateVcf.par_bed": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed",
  "AnnotateVcf.prefix": "dfci-ufc.v1.sv.annotation_test.chr1",
  "AnnotateVcf.protein_coding_gtf": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/gencode.v47.basic.protein_coding.canonical.gtf",
  "AnnotateVcf.runtime_attr_svannotate": {"cpu_cores" : 4, "mem_gb": 7.5},
  "AnnotateVcf.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "AnnotateVcf.sv_per_shard": 5000,
  "AnnotateVcf.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-09-02-v1.0.5-631368eb",
  "AnnotateVcf.vcf": "$WORKSPACE_BUCKET/scratch/dfci-g2c.v1.sv.redundant_collapsed.cpx_only.chr1.vcf.gz"
}
EOF

# Submit chr1 annotation test
cromshell --no_turtle -t 120 -mc submit \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip svannotation.dependencies.zip \
  code/wdl/gatk-sv/svannotation_v1.1/AnnotateVcf.wdl \
  cromshell/inputs/AnnotateVcf.inputs.chr1_test.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/annotation_test.job_ids.list

# Monitor annotation test
monitor_workflow $( tail -n1 cromshell/job_ids/annotation_test.job_ids.list )


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


#############################################
# Collect GATK-SV QC prior to GT imputation #
#############################################

# Note: this workflow below is scattered across all five workspaces for 
# max parallelization. It must be submitted as below in each workspace.

# Reaffirm staging directory
staging_dir=staging/initial_gatksv_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectInitialGatksvQcMetrics.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 250,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.001,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:1aac84d",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 5000,
  "CollectVcfQcMetrics.output_prefix": "dfci-g2c.v1.initial_gatksv_qc.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.sample_priority_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.sample_qc_priority.tsv",
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/\$CONTIG/HardFilterPart2/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.identical.reclustered.posthoc_filtered.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/\$CONTIG/HardFilterPart2/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.identical.reclustered.posthoc_filtered.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metric collection workflows
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectInitialGatksvQcMetrics.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/GatksvQcMetrics/ \
  --name CollectInitialGatksvQcMetrics \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.CollectInitialGatksvQcMetrics.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --submission-gate 3 \
  --max-attempts 3


#############################
# Curate GATK-SV QC targets #
#############################

# Reaffirm staging directory
staging_dir=staging/qc_targets
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Estimate number of variants per genome in gnomAD for the necessary contigs
for k in $( seq 1 22 ) X Y; do
  gsutil cat \
    gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/chr$k/gnomad.v4.1.gatksv.chr$k.*.sites.bed.gz
done | gunzip -c \
| code/scripts/estimate_vpg_from_sites.py \
| fgrep -v "#" \
| awk -v OFS="\t" '{ print "variants_per_genome."$1":median", $2 }' \
> $staging_dir/dfci-g2c.v1.gatksv.qc_targets.tsv

# Copy QC targets to central bucket for reference by Cromwell
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.gatksv.qc_targets.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv-gt-imputation-qc/


##################################################
# Analyze & visualize initial GATK-SV QC metrics #
##################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/initial_gatksv_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

cat << EOF > $staging_dir/main_keys.list
size_distrib
af_distrib
size_vs_af_distrib
all_svs_bed
common_svs_bed
genotype_distrib
ld_stats
EOF

cat << EOF > $staging_dir/bench_keys.list
site_benchmark_ppv_by_freqs
site_benchmark_sensitivity_by_freqs
site_benchmark_common_sv_ppv_beds
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
  json_fname=CollectInitialGatksvQcMetrics.chr$k.outputs.json
  gsutil cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/GatksvQcMetrics/chr$k/$json_fname \
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
    echo "gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/chr$k/gnomad.v4.1.gatksv.chr$k.$suffix.merged.tsv.gz" \
    >> $staging_dir/gnomAD_$suffix.uris.list
  done
done

# Write input .json
cat << EOF | python -m json.tool > cromshell/inputs/PlotInitialGatksvQcMetrics.inputs.json
{
  "PlotVcfQcMetrics.af_distribution_tsvs": $( collapse_txt $staging_dir/af_distrib.uris.list ),
  "PlotVcfQcMetrics.all_sv_beds": $( collapse_txt $staging_dir/all_svs_bed.uris.list ),
  "PlotVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PlotVcfQcMetrics.benchmark_interval_names": ["Easy", "Hard"],
  "PlotVcfQcMetrics.common_af_cutoff": 0.001,
  "PlotVcfQcMetrics.common_sv_beds": $( collapse_txt $staging_dir/common_svs_bed.uris.list ),
  "PlotVcfQcMetrics.custom_qc_target_metrics": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv-gt-imputation-qc/dfci-g2c.v1.gatksv.qc_targets.tsv",
  "PlotVcfQcMetrics.deduplicate": true,
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:cd2ca89",
  "PlotVcfQcMetrics.output_prefix": "dfci-g2c.v1.initial_gatksv_qc",
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
  "PlotVcfQcMetrics.site_benchmark_common_sv_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_easy.uris.list ),
                                                           $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_hard.uris.list ) ]],
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
  cromshell/inputs/PlotInitialGatksvQcMetrics.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.PlotInitialGatksvQcMetrics.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotInitialGatksvQcMetrics.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotGatksvQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotInitialGatksvQcMetrics.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotGatksvQc/

# Clear Cromwell execution & output buckets for plotting job
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.PlotInitialGatksvQcMetrics.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwel*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


#################################################
# Prepare covariates for SV genotype imputation #
#################################################

# This block only needs to be run once for the entire cohort

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


#################################################
# Refine common SV genotypes with flanking SNVs #
#################################################

# This block must be run once per workspace for parallelization

# Reaffirm staging directory
staging_dir=staging/sv_gt_cleanup
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Make lists of all SNV VCFs and indexes for each contig
while read contig; do
  gsutil -m ls \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart2/$contig/**vcf.gz \
  | sort -V \
  | awk -v OFS="\t" '{ print $1, $1".tbi" }' \
  > $staging_dir/dfci-g2c.v1.sv_regenotyping.snv_vcf_info.$contig.tsv
done < contig_lists/dfci-g2c.v1.contigs.$WN.list
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.sv_regenotyping.snv_vcf_info.*.tsv \
  $MAIN_WORKSPACE_BUCKET/data/sv_regenotyping/

# Write template input .json for SV GT refinement
cat << EOF > $staging_dir/RefineSvGenotypesWithSnvs.inputs.template.json
{
  "RefineSvGenotypesWithSnvs.ConcatVcfs.boot_disk_gb": 25,
  "RefineSvGenotypesWithSnvs.ConcatVcfs.disk_gb": 500,
  "RefineSvGenotypesWithSnvs.ConcatVcfs.mem_gb": 7.5,
  "RefineSvGenotypesWithSnvs.QuerySnvs.n_preemptible": 0,
  "RefineSvGenotypesWithSnvs.g2c_analysis_docker": "vanallenlab/g2c_analysis:a9d85cd",
  "RefineSvGenotypesWithSnvs.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "RefineSvGenotypesWithSnvs.linux_docker": "ubuntu:plucky-20251001",
  "RefineSvGenotypesWithSnvs.min_imputation_r2": 0.2,
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
  "RefineSvGenotypesWithSnvs.snv_vcfs_per_shard": 125,
  "RefineSvGenotypesWithSnvs.svs_per_shard": 125,
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
  --submission-gate 240 \
  --vm-gate 600 \
  --max-attempts 3


##############################################
# Fix typo in SV VCF header after imputation #
##############################################

# This is only necessary because there was a typo in FORMAT/IMP
# This should not be necessary for all future cohorts (typo/bug was fixed)

# This is a lightweight WDL so can easily be parallelized for all chromosomes in a single workspace

# Reaffirm staging directory
staging_dir=staging/sv_gt_cleanup
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json
cat << EOF > $staging_dir/FixHeaderTypo.inputs.template.json
{
  "FixHeaderTypo.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "FixHeaderTypo.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.imputed.vcf.gz",
  "FixHeaderTypo.vcf_idx": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.imputed.vcf.gz.tbi"
}
EOF

# Submit, monitor, and stage/cleanup
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/FixHeaderTypo.wdl \
  --input-json-template $staging_dir/FixHeaderTypo.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup_header_fix \
  --status-tsv cromshell/progress/dfci-g2c.v1.FixHeaderTypo.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 10 \
  --submission-gate 0 \
  --max-attempts 2


##########################################
# Collect GATK-SV QC after GT imputation #
##########################################

# Note: this workflow below is scattered across all five workspaces for 
# max parallelization. It must be submitted as below in each workspace.

# Reaffirm staging directory
staging_dir=staging/gatksv_qc_post_imputation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectGatksvQcPostImputation.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 250,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.001,
  "CollectVcfQcMetrics.extra_vcf_preprocessing_commands": "| bcftools view -i 'AC > 0 | FILTER = \"MULTIALLELIC\"'",
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:a9d85cd",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 5000,
  "CollectVcfQcMetrics.output_prefix": "dfci-g2c.v1.gatksv_qc_post_imputation.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.sample_priority_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.sample_qc_priority.tsv",
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup_header_fix/\$CONTIG/FixTypo/dfci-g2c.v1.\$CONTIG.imputed.typo_fixed.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup_header_fix/\$CONTIG/FixTypo/dfci-g2c.v1.\$CONTIG.imputed.typo_fixed.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metric collection workflows
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectGatksvQcPostImputation.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv-gt-imputation-qc/GatksvQcMetrics/ \
  --name CollectGatksvQcPostImputation \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.CollectGatksvQcPostImputation.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --submission-gate 3 \
  --max-attempts 3


##############################################################
# Analyze & visualize GATK-SV QC metrics pre/post imputation #
##############################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/gatksv_qc_post_imputation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

cat << EOF > $staging_dir/main_keys.list
size_distrib
af_distrib
size_vs_af_distrib
all_svs_bed
common_svs_bed
genotype_distrib
ld_stats
EOF

cat << EOF > $staging_dir/bench_keys.list
site_benchmark_ppv_by_freqs
site_benchmark_sensitivity_by_freqs
site_benchmark_common_sv_ppv_beds
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
  json_fname=CollectGatksvQcPostImputation.chr$k.outputs.json
  gsutil cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv-gt-imputation-qc/GatksvQcMetrics/chr$k/$json_fname \
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
    echo "gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/chr$k/gnomad.v4.1.gatksv.chr$k.$suffix.merged.tsv.gz" \
    >> $staging_dir/gnomAD_$suffix.uris.list
  done
done

# Write input .json
cat << EOF | python -m json.tool > cromshell/inputs/PlotGatksvQcPostImputation.inputs.json
{
  "PlotVcfQcMetrics.af_distribution_tsvs": $( collapse_txt $staging_dir/af_distrib.uris.list ),
  "PlotVcfQcMetrics.all_sv_beds": $( collapse_txt $staging_dir/all_svs_bed.uris.list ),
  "PlotVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PlotVcfQcMetrics.benchmark_interval_names": ["Easy", "Hard"],
  "PlotVcfQcMetrics.common_af_cutoff": 0.001,
  "PlotVcfQcMetrics.common_sv_beds": $( collapse_txt $staging_dir/common_svs_bed.uris.list ),
  "PlotVcfQcMetrics.custom_qc_target_metrics": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv-gt-imputation-qc/dfci-g2c.v1.gatksv.qc_targets.tsv",
  "PlotVcfQcMetrics.deduplicate": true,
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:cd2ca89",
  "PlotVcfQcMetrics.output_prefix": "dfci-g2c.v1.gatksv_qc_post_imputation",
  "PlotVcfQcMetrics.peak_ld_stat_tsvs": $( collapse_txt $staging_dir/ld_stats.uris.list ),
  "PlotVcfQcMetrics.PlotSiteBenchmarking.mem_gb": 32,
  "PlotVcfQcMetrics.PlotSiteBenchmarking.n_cpu": 8,
  "PlotVcfQcMetrics.PlotSiteMetrics.mem_gb": 32,
  "PlotVcfQcMetrics.PlotSiteMetrics.n_cpu": 8,
  "PlotVcfQcMetrics.previous_stats": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotGatksvQc/dfci-g2c.v1.initial_gatksv_qc.all_qc_summary_metrics.tsv",
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
  "PlotVcfQcMetrics.site_benchmark_common_sv_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_easy.uris.list ),
                                                           $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_hard.uris.list ) ]],
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
  cromshell/inputs/PlotGatksvQcPostImputation.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.PlotGatksvQcPostImputation.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotGatksvQcPostImputation.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv-gt-imputation-qc/PlotGatksvQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotGatksvQcPostImputation.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv-gt-imputation-qc/PlotGatksvQc/

# Clear Cromwell execution & output buckets for plotting job
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.PlotGatksvQcPostImputation.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwel*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


####################################
# Determine final variant sharding #
####################################

# The below must be run once for each workspace

# Reaffirm staging directory
staging_dir=staging/indel_sv_integration
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi
gsutil cp gs://dfci-g2c-refs/hg38/hg38.genome $staging_dir/

# Download & index raw VCF QC maps
while read contig; do
  # Prep contig-specific directory
  csdir=$staging_dir/$contig
  if [ -e $csdir ]; then rm -rf $csdir; fi
  mkdir $csdir

  # Localize & index variant maps
  for vc in snvs indels svs; do
    key="all_${vc}_bed"
    gsutil -m cat \
      $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/VcfQcMetrics/$contig/CollectInitialVcfQcMetrics.$contig.outputs.json \
    | jq .\"CollectVcfQcMetrics.$key\" \
    | tr -d '"' \
    | gsutil -m cp -I $csdir/
    tabix -p bed -f $csdir/dfci-g2c.v1.initial_qc.$contig.all_$vc.bed.gz
  done
done < contig_lists/dfci-g2c.v1.contigs.$WN.list

# Define final analysis shard intervals for SNVs, indels, and SVs
# Rough logic: all chromosomes in a workspace should sum to ~2x AoU Cromwell quota (1.1k x 2 ~ 2k)
# These 2k shards sholud be divided among chromosomes based on total variant count
# Then, within each chromosome, they should be partitioned by variant class count
# And along each chromosome should be segmented based on density
# Desired end result: all VCF shards should have roughly the same number of records

# First, get variant counts by class per contig
while read contig; do
  for wrapper in 1; do
    echo $contig
    for vc in snv indel sv; do
      zcat $staging_dir/$contig/dfci-g2c.v1.initial_qc.$contig.all_${vc}s.bed.gz \
      | grep -ve '^#' | wc -l
    done
  done | paste -s \
  | awk -v FS="\t" -v OFS="\t" '{ sum=$2+$3+$4 }END{ print $0, sum }'
done < contig_lists/dfci-g2c.v1.contigs.$WN.list \
> $staging_dir/contig.variant_counts.tsv

# Second, determine shards allocated per contig
denom=$( awk '{ sum+=$5 }END{ print sum }' $staging_dir/contig.variant_counts.tsv )
awk -v scalar=2000 -v denom=$denom -v FS="\t" -v OFS="\t" \
  '{ print $1, int(scalar * $5 / denom) }' \
  $staging_dir/contig.variant_counts.tsv \
> $staging_dir/shards_per_contig.tsv

# Third, partition shards across variant classes per contig
while read contig; do
  # Compute number of variants to allocate per shard
  total_var=$( awk -v FS="\t" -v contig=$contig \
                 '{ if ($1==contig) print $5 }' \
                 $staging_dir/contig.variant_counts.tsv )
  total_shards=$( awk -v FS="\t" -v contig=$contig \
                    '{ if ($1==contig) print $2 }' \
                    $staging_dir/shards_per_contig.tsv )
  vps=$( echo "" | awk -v n=$total_var -v d=$total_shards '{ print int(n/d) }' )
  
  # Shard intervals for each variant class
  csdir=$staging_dir/$contig
  awk -v contig=$contig -v OFS="\t" \
    '{ if ($1==contig) print contig, 1, $2, "+", contig }' \
    $staging_dir/hg38.genome \
  > $csdir/$contig.full.interval_list
  for vc in snv indel sv; do
    code/scripts/split_intervals.py \
      -i $csdir/$contig.full.interval_list \
      --var-sites $csdir/dfci-g2c.v1.initial_qc.$contig.all_${vc}s.bed.gz \
      --vars-per-shard $vps \
      --bed-style \
      --verbose \
    | awk -v OFS="\t" -v prefix="dfci-g2c.v1.$vc.$contig" \
      '{ print $0, prefix"."NR }' \
    | bgzip -c \
    > $staging_dir/dfci-g2c.v1.analysis_shards.$contig.$vc.bed.gz
    tabix -p bed -f $staging_dir/dfci-g2c.v1.analysis_shards.$contig.$vc.bed.gz
  done
done < contig_lists/dfci-g2c.v1.contigs.$WN.list

# Once complete, copy all final sharded intervals to a permanent staging bucket
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.analysis_shards.chr*.*.bed.gz* \
  $MAIN_WORKSPACE_BUCKET/data/g2c_partition_maps/


########################################
# Integrate small SVs and large indels #
########################################

# Reaffirm staging directory
staging_dir=staging/indel_sv_integration
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Curate reference data required for this workflow
# This only need to be run once for the project
# This should be run locally so it can be staged in a public bucket
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
zcat gap.txt.gz \
| cut -f2-4 \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| fgrep -v "_" \
| grep -e '^chr' \
| bgzip -c \
> hg38.gaps.bed.gz
gsutil -m cp hg38.gaps.bed.gz gs://dfci-g2c-refs/hg38/
gsutil -m cp gs://dfci-g2c-refs/hg38/hg38.genome ./
mkdir contig_genome_files
for k in $( seq 1 22 ) X Y; do
  contig="chr$k"
  awk -v contig=$contig '{ if ($1==contig) print }' hg38.genome \
  > contig_genome_files/hg38.$contig.genome
done
gsutil -m cp -r contig_genome_files gs://dfci-g2c-refs/hg38/

# All of the below must be run once for each workspace

# Write template input .json 
cat << EOF > $staging_dir/UnifyGatkCallsets.inputs.template.json
{
  "UnifyGatkCallsets.g2c_analysis_docker": "vanallenlab/g2c_analysis:11983b2",
  "UnifyGatkCallsets.g2c_analysis_docker_dev_tmp": "vanallenlab/g2c_analysis:691fc72",
  "UnifyGatkCallsets.gatkhc_vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/data/sv_regenotyping/dfci-g2c.v1.sv_regenotyping.snv_vcf_info.\$CONTIG.tsv",
  "UnifyGatkCallsets.gatksv_vcfs": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup_header_fix/\$CONTIG/FixTypo/dfci-g2c.v1.\$CONTIG.imputed.typo_fixed.vcf.gz"],
  "UnifyGatkCallsets.gatksv_vcf_idxs": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup_header_fix/\$CONTIG/FixTypo/dfci-g2c.v1.\$CONTIG.imputed.typo_fixed.vcf.gz.tbi"],
  "UnifyGatkCallsets.genome_file": "gs://dfci-g2c-refs/hg38/contig_genome_files/hg38.\$CONTIG.genome",
  "UnifyGatkCallsets.indel_partition_intervals": "$MAIN_WORKSPACE_BUCKET/data/g2c_partition_maps/dfci-g2c.v1.analysis_shards.\$CONTIG.indel.bed.gz",
  "UnifyGatkCallsets.nmask_bed": "gs://dfci-g2c-refs/hg38/hg38.gaps.bed.gz",
  "UnifyGatkCallsets.PartitionSnvOutputs.n_preemptible": 0,
  "UnifyGatkCallsets.snv_partition_intervals": "$MAIN_WORKSPACE_BUCKET/data/g2c_partition_maps/dfci-g2c.v1.analysis_shards.\$CONTIG.snv.bed.gz",
  "UnifyGatkCallsets.SplitIndelsForClustering.mem_gb": 32,
  "UnifyGatkCallsets.SplitIndelsForClustering.n_cpu": 16,
  "UnifyGatkCallsets.SplitIndelsForClustering.n_preemptible": 0,
  "UnifyGatkCallsets.sv_partition_intervals": "$MAIN_WORKSPACE_BUCKET/data/g2c_partition_maps/dfci-g2c.v1.analysis_shards.\$CONTIG.sv.bed.gz"
}
EOF

# DEV NOTE: ONLY CHR22 FOR NOW

# Submit, monitor, and stage/cleanup SV GT refinement
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/UnifyGatkCallsets.wdl \
  --input-json-template $staging_dir/UnifyGatkCallsets.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/indel_sv_integration \
  --contig-list <( echo "chr22" ) \
  --status-tsv cromshell/progress/dfci-g2c.v1.UnifyGatkCallsets.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 10 \
  --submission-gate 0.1 \
  --max-attempts 3 \
  --dry-run --hard-reset

# DEV: submit
cromshell -t 120 abort $( tail -n1 scratch/dev.wid.list )
cromshell --no_turtle -t 120 -mc submit \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip g2c.dependencies.zip --no-validation \
  code/wdl/pancan_germline_wgs/UnifyGatkCallsets.wdl \
  /home/jupyter/cromshell/inputs/UnifyGatkCallsets.inputs.chr22.json \
| jq .id | tr -d '"' \
>> scratch/dev.wid.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 scratch/dev.wid.list ) 2



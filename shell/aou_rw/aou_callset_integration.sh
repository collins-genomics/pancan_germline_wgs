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
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.snv.sites.bed.gz"],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.indel.sites.bed.gz"],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.sv.sites.bed.gz"],
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
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:1aac84d",
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
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotQc/

# Clear Cromwell execution & output buckets for plotting job
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.PlotInitialGatksvQcMetrics.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell/*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


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
  "RefineSvGenotypesWithSnvs.QuerySnvs.n_preemptible": 0,
  "RefineSvGenotypesWithSnvs.g2c_analysis_docker": "vanallenlab/g2c_analysis:1aac84d",
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
  --max-attempts 3


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
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:1aac84d",
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
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.snv.sites.bed.gz"],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.indel.sites.bed.gz"],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.imputed.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.imputed.vcf.gz.tbi"]
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
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:1aac84d",
  "PlotVcfQcMetrics.output_prefix": "dfci-g2c.v1.gatksv_qc_post_imputation",
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
  cromshell/inputs/PlotGatksvQcPostImputation.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.PlotGatksvQcPostImputation.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotGatksvQcPostImputation.job_ids.list ) 5

# # Once workflow is complete, stage output
# gsutil -m rm -rf $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotQc
# cromshell -t 120 list-outputs \
#   $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotGatksvQcPostImputation.job_ids.list ) \
# | awk '{ print $2 }' \
# | gsutil -m cp -I \
#   $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotQc/

# # Clear Cromwell execution & output buckets for plotting job
# gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.PlotGatksvQcPostImputation.job_ids.list \
#                 | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell/*/PlotVcfQcMetrics/" \
#                   '{ print bucket_prefix$1"/**" }' ) \
# > uris_to_delete.list
# cleanup_garbage


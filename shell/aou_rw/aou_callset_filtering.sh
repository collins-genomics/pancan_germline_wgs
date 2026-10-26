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


####################################################################
# [Alpha] Apply pre-trained high-sensitivity FilterGenotypes model #
####################################################################

# Make staging directory
staging_dir=staging/OnlyFilterGenotypesParallel
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write text file of SV VCF URIs for each contig
while read contig; do
  gsutil -m ls \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/indel_sv_integration/$contig/**sv*.vcf.gz \
  | sort -Vk1,1 \
  > $staging_dir/dfci-g2c.v1.filtergenotypes_qc.sv_vcf_info.$contig.txt
done < contig_lists/dfci-g2c.v1.contigs.$WN.list

# Initialize .json of contig-specific VCF array overrides
echo "{ " > $staging_dir/OnlyFilterGenotypesParallel.contig_variable_overrides.json
while read contig; do
  echo "\"$contig\" : {\"VCF_ARRAY\" : $( collapse_txt $staging_dir/dfci-g2c.v1.filtergenotypes_qc.sv_vcf_info.$contig.txt )},"
done < contig_lists/dfci-g2c.v1.contigs.$WN.list \
| paste -s -d\  | sed 's/,$//g' \
>> $staging_dir/OnlyFilterGenotypesParallel.contig_variable_overrides.json
echo " }" >> $staging_dir/OnlyFilterGenotypesParallel.contig_variable_overrides.json

# Write template .json for module 19
cat << EOF > $staging_dir/OnlyFilterGenotypesParallel.inputs.template.json
{
  "OnlyFilterGenotypesParallel.filter_script": "$MAIN_WORKSPACE_BUCKET/code/scripts/apply_sl_filter.g2c_imputed_integrated.py",
  "OnlyFilterGenotypesParallel.FilterVcfs.filter_machine_type": "n2d-standard-4",
  "OnlyFilterGenotypesParallel.FilterVcfs.sanitize_machine_type": "n2d-standard-4",
  "OnlyFilterGenotypesParallel.FilterVcfs.shard_vcf": false,
  "OnlyFilterGenotypesParallel.no_call_rate_cutoff": 1,
  "OnlyFilterGenotypesParallel.ploidy_table": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/17/dfci-g2c.v1.ploidy.tsv",
  "OnlyFilterGenotypesParallel.sl_filter_args": "--small-del-threshold -53 --medium-del-threshold -12 --small-dup-threshold -105 --medium-dup-threshold -81 --ins-threshold -97",
  "OnlyFilterGenotypesParallel.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "OnlyFilterGenotypesParallel.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-06-27-v1.0.4-63e6c81e",
  "OnlyFilterGenotypesParallel.vcfs": \$VCF_ARRAY
}
EOF

# Run filtering component of GATK-SV module 19 on all SV VCFs
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/OnlyFilterGenotypesParallel.wdl \
  --input-json-template $staging_dir/OnlyFilterGenotypesParallel.inputs.template.json \
  --contig-variable-overrides $staging_dir/OnlyFilterGenotypesParallel.contig_variable_overrides.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/data/alpha_v0.9/OnlyFilterGenotypesParallel \
  --dependencies-zip gatksv.dependencies.zip \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/OnlyFilterGenotypesParallel.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v0.9" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 5


#####################################################################
# [Alpha] Collect QC metrics after high-sensitivity FilterGenotypes #
#####################################################################

# Note: this workflow below is scattered across all five workspaces for 
# max parallelization. It must be submitted as below in each workspace.

# Rotate Cromwell cache before embarking on these workflows, which have large scatter counts
~/code/scripts/rotate_cromwell_cache.sh delete

# Reaffirm staging directory
staging_dir=staging/filtergenotypes_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write two-column .tsv of VCF & index info for each contig
while read contig; do
  gsutil cat \
    $MAIN_WORKSPACE_BUCKET/data/alpha_v0.9/OnlyFilterGenotypesParallel/$contig/OnlyFilterGenotypesParallel.$contig.outputs.json \
  | jq '.["OnlyFilterGenotypesParallel.filtered_vcfs"]' \
  | fgrep "gs://" | awk '{ print $1 }' | tr -d '",' \
  | awk -v OFS="\t" '{ print $1, $1".tbi" }' \
  > $staging_dir/dfci-g2c.v0.9.post_filtergenotypes_qc.vcf_info.$contig.tsv
done < contig_lists/dfci-g2c.v1.contigs.$WN.list
gsutil -m cp \
  $staging_dir/dfci-g2c.v0.9.post_filtergenotypes_qc.vcf_info.*.tsv \
  $MAIN_WORKSPACE_BUCKET/data/alpha_v0.9/filtergenotypes_qc/vcf_list_inputs/

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectFilterGenotypesQc.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 250,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.001,
  "CollectVcfQcMetrics.extra_vcf_preprocessing_commands": "| bcftools view -i 'AC > 0 | FILTER = \"MULTIALLELIC\"'",
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 5000,
  "CollectVcfQcMetrics.output_prefix": "dfci-g2c.v0.9.sv_qc_post_filtergenotypes.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.gcp_machine_type": "n2d-standard-4",
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
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/data/alpha_v0.9/filtergenotypes_qc/vcf_list_inputs/dfci-g2c.v0.9.post_filtergenotypes_qc.vcf_info.\$CONTIG.tsv"
}
EOF

# Submit, monitor, stage, and cleanup QC metric collection workflows
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectFilterGenotypesQc.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/data/alpha_v0.9/CollectFilterGenotypesQc/ \
  --name CollectFilterGenotypesQc \
  --contig-list contig_lists/dfci-g2c.v1.contigs.$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v0.9.CollectFilterGenotypesQc.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v0.9" \
  --outer-gate 30 \
  --submission-gate 240 \
  --max-attempts 3

# Clear Cromwell cache after finishing these workflows
~/code/scripts/rotate_cromwell_cache.sh delete


#######################################################################
# [Alpha] Visualize QC metrics after high-sensitivity FilterGenotypes #
#######################################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/filtergenotypes_qc
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
  json_fname=CollectFilterGenotypesQc.chr$k.outputs.json
  gsutil cp \
    $MAIN_WORKSPACE_BUCKET/data/alpha_v0.9/CollectFilterGenotypesQc/chr$k/$json_fname \
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
cat << EOF | python -m json.tool > cromshell/inputs/PlotFilterGenotypesQc.inputs.json
{
  "PlotVcfQcMetrics.af_distribution_tsvs": $( collapse_txt $staging_dir/af_distrib.uris.list ),
  "PlotVcfQcMetrics.all_sv_beds": $( collapse_txt $staging_dir/all_svs_bed.uris.list ),
  "PlotVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PlotVcfQcMetrics.benchmark_interval_names": ["Easy", "Hard"],
  "PlotVcfQcMetrics.common_af_cutoff": 0.001,
  "PlotVcfQcMetrics.common_sv_beds": $( collapse_txt $staging_dir/common_svs_bed.uris.list ),
  "PlotVcfQcMetrics.custom_qc_target_metrics": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.qc_targets.tsv",
  "PlotVcfQcMetrics.deduplicate": true,
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "PlotVcfQcMetrics.output_prefix": "dfci-g2c.v0.9.filtergenotypes_qc",
  "PlotVcfQcMetrics.peak_ld_stat_tsvs": $( collapse_txt $staging_dir/ld_stats.uris.list ),
  "PlotVcfQcMetrics.PlotSiteBenchmarking.gcp_machine_type": "n2d-standard-8",
  "PlotVcfQcMetrics.PlotSiteMetrics.gcp_machine_type": "n2d-standard-16",
  "PlotVcfQcMetrics.previous_stats": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/integrated-indel-sv-qc/PlotDevQc/dfci-g2c.v1.integrated_qc.all_qc_summary_metrics.tsv",
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
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotFilterGenotypesQc.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v0.9.PlotFilterGenotypesQc.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v0.9.PlotFilterGenotypesQc.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $MAIN_WORKSPACE_BUCKET/data/alpha_v0.9/PlotFilterGenotypesQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v0.9.PlotFilterGenotypesQc.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/data/alpha_v0.9/PlotFilterGenotypesQc/

# Clear Cromwell execution & output buckets for plotting job
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v0.9.PlotFilterGenotypesQc.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/workflows/cromwell*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


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



#########################################################
# [Main] Curate sample-level features for all filtering #
#########################################################

# TODO: implement this


# ###################################################################
# # Curate variant- and genotype-level features for SV GT filtering #
# ###################################################################

# # BELOW IS DEV: manual submit for chr22 to familiarize with wb workflow syntax

# # Reaffirm staging directory
# staging_dir=~/staging/sv_gt_feature_curation
# if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# # Dev: set contig manually (this can be looped once automated)
# contig=chr22

# # Create VCF input array
# # NOTE: because this was the first task after migrating from AoU RW v1.0 to Verily,
# # we need to update the URI prefixes for all VCFs to the new Verily buckets
# gsutil -m cat \
#   $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/indel_sv_integration/$contig/UnifyGatkCallsets.$contig.outputs.json \
# | jq .\"UnifyGatkCallsets.cleaned_sv_vcfs\" \
# | fgrep "gs://" | tr -d '",' \
# | awk -v FS="42de3578da45" -v prefix="$MAIN_WORKSPACE_BUCKET" -v OFS="\t" \
#   '{ print prefix$2, prefix$2".tbi" }' \
# > $staging_dir/CollectGTFilterFeatures.vcf_inputs.$contig.tsv
# gsutil -m cp \
#   $staging_dir/CollectGTFilterFeatures.vcf_inputs.$contig.tsv \
#   $MAIN_WORKSPACE_BUCKET/data/sv_gt_filtering/

# # # Create inputs and upload to workspace bucket
# # input_json="CollectGTFilterFeatures.$contig.inputs.json"
# # cat << EOF > "cromwell/inputs/$input_json"
# # {
# #   "CollectGTFilterFeatures.bed_features": ["gs://dfci-g2c-refs/giab/$contig/giab.hg38.broad_callable.hard.$contig.bed.gz",
# #                                            "gs://dfci-g2c-refs/ucsc/hg38/hg38.segdups.bed.gz",
# #                                            "gs://dfci-g2c-refs/ucsc/hg38/hg38.simple_repeats.bed.gz"],
# #   "CollectGTFilterFeatures.bigwig_features": ["gs://dfci-g2c-refs/ucsc/hg38/hg38.100mer_multiUmap.bw"],
# #   "CollectGTFilterFeatures.bigwig_feature_names": ["umap"],
# #   "CollectGTFilterFeatures.bed_feature_names": ["giab_hard", "segdup", "simrep"],
# #   "CollectGTFilterFeatures.g2c_analysis_docker": "vanallenlab/g2c_analysis:77bb125",
# #   "CollectGTFilterFeatures.linux_docker": "ubuntu:plucky-20251001",
# #   "CollectGTFilterFeatures.output_prefix": "dfci-g2c.v1.$contig.sv",
# #   "CollectGTFilterFeatures.vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/data/sv_gt_filtering/CollectGTFilterFeatures.vcf_inputs.$contig.tsv"
# # }
# # EOF
# # gsutil cp "cromwell/inputs/$input_json" $WORKSPACE_BUCKET/cromwell-inputs/$workflow_base/

# # # Test submission for a single chr22 SV VCF
# # wb workflow job run \
# #   --workflow "CollectGTFilterFeatures" \
# #   --inputs-uri "$WORKSPACE_BUCKET/cromwell-inputs/$workflow_base/$input_json" \
# #   --output-bucket-id "$( echo $WORKSPACE_BUCKET | sed 's/^gs\:\/\///g'  )" \
# #   --output-path cromwell-execution \
# #   --read-from-cache \
# #   --write-to-cache \
# #   --format JSON \
# # | python -m json.tool \
# # > "cromwell/submissions/CollectGTFilterFeatures.$contig.submission.json"

# # Write template input .json for SV GT refinement
# cat << EOF > $staging_dir/CollectSVGTFilterFeatures.inputs.template.json
# {
#   "CollectGTFilterFeatures.bed_features": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz",
#                                            "gs://dfci-g2c-refs/ucsc/hg38/hg38.segdups.bed.gz",
#                                            "gs://dfci-g2c-refs/ucsc/hg38/hg38.simple_repeats.bed.gz"],
#   "CollectGTFilterFeatures.bigwig_features": ["gs://dfci-g2c-refs/ucsc/hg38/hg38.100mer_multiUmap.bw"],
#   "CollectGTFilterFeatures.bigwig_feature_names": ["umap"],
#   "CollectGTFilterFeatures.bed_feature_names": ["giab_hard", "segdup", "simrep"],
#   "CollectGTFilterFeatures.g2c_analysis_docker": "vanallenlab/g2c_analysis:77bb125",
#   "CollectGTFilterFeatures.linux_docker": "ubuntu:plucky-20251001",
#   "CollectGTFilterFeatures.output_prefix": "dfci-g2c.v1.\$CONTIG.sv",
#   "CollectGTFilterFeatures.vcf_info_tsv": "$MAIN_WORKSPACE_BUCKET/data/sv_gt_filtering/CollectGTFilterFeatures.vcf_inputs.\$CONTIG.tsv"
# }
# EOF

# # Submit, monitor, and stage/cleanup SV GT refinement
# code/scripts/manage_chromshards.py \
#   --wdl code/wdl/pancan_germline_wgs/CollectGTFilterFeatures.wdl \
#   --input-json-template $staging_dir/CollectSVGTFilterFeatures.inputs.template.json \
#   --dependencies-zip g2c.dependencies.zip \
#   --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/CollectSVGTFilterFeatures \
#   --name CollectSVGTFilterFeatures \
#   --contig-list <( fgrep chr22 contig_lists/dfci-g2c.v1.contigs.$WN.list ) \
#   --status-tsv cromshell/progress/dfci-g2c.v1.CollectSVGTFilterFeatures.progress.tsv \
#   --workflow-id-log-prefix "dfci-g2c.v1"



# # # Check progress
# # wid=$( jq .runId "cromwell/submissions/CollectGTFilterFeatures.$contig.submission.json" | tr -d '"' )
# # monitor_workflow $wid
# # wb workflow job describe --job-id=$wid --format=JSON | python -m json.tool

# # TODO: finish implementing this



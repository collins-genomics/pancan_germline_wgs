#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to perform intake QC & batching for DFCI G2C phase 1 cohort

# Note that this code is designed to be run inside the AoU Researcher Workbench


# Note2: this code was defined for the AoU RW v1.0, and it may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in data data/tarballs data/gatksv_batch_info data/gatksv_batch_info/sample_lists \
           plots plots/raw_intake_qc plots/global_qc_pass plots/batch_qc_pass; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
for suffix in py R sh; do
  find code/ -name "*.${suffix}" | xargs -I {} chmod a+x {}
done

# Source .bashrc
. ~/code/refs/dotfiles/aou.rw.bashrc

# Install necessary packages
. code/refs/install_packages.sh R

# Localize intake QC data and references
gsutil -m -u $GPROJECT cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/dfci-g2c.intake_qc.aou.tsv.gz \
  gs://dfci-g2c-inputs/intake_qc/dfci-g2c.intake_qc.non_aou.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/phenotypes/dfci-g2c.aou.phenos.tsv.gz \
  gs://dfci-g2c-inputs/phenotypes/dfci-g2c.non_aou.phenos.tsv.gz \
  gs://dfci-g2c-refs/hgsv/HGSV.ByrskaBishop.sample_populations.tsv \
  data/
gsutil -m -u $GPROJECT cat \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/other/all_g2c.phase1.samples.list \
  gs://dfci-g2c-inputs/intake_qc/all_g2c.phase1.samples.non_aou.list \
| sort -V > data/all_g2c.phase1.intake_qc.samples.list


########################
# INTAKE QC & BATCHING #
########################

# Add minimal necessary phenotype data to non-AoU and AoU samples for QC + batching
for subset in non_aou aou; do
  code/scripts/join_qc_and_phenotypes.R \
    --qc-tsv data/dfci-g2c.intake_qc.$subset.tsv.gz \
    --phenos-tsv data/dfci-g2c.$subset.phenos.tsv.gz \
    --out-tsv data/dfci-g2c.intake_qc_with_phenos.$subset.tsv
  gzip -f data/dfci-g2c.intake_qc_with_phenos.$subset.tsv
done

# Merge AoU and non-AoU manifests
# Also labels ancestry, sex, and infers PCR status
code/scripts/merge_intake_qc_manifests.R \
  --aou-tsv data/dfci-g2c.intake_qc_with_phenos.aou.tsv.gz \
  --non-aou-tsv data/dfci-g2c.intake_qc_with_phenos.non_aou.tsv.gz \
  --include-samples data/all_g2c.phase1.intake_qc.samples.list \
  --hgsv-pop-assignments data/HGSV.ByrskaBishop.sample_populations.tsv \
  --id-prefix G2C \
  --suffix-length 6 \
  --out-tsv data/dfci-g2c.intake_qc.all.tsv \
  --out-ped data/dfci-g2c.all_samples.ped
gzip -f data/dfci-g2c.intake_qc.all.tsv

# Copy merged manifest to main workspace bucket
gsutil -m cp \
  data/dfci-g2c.intake_qc.all.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/

# Copy .ped file to GATK-SV bucket
gsutil -m cp \
  data/dfci-g2c.all_samples.ped \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/

# Generate initial QC plots prior to batching
code/scripts/plot_intake_qc.R \
  --qc-tsv data/dfci-g2c.intake_qc.all.tsv.gz \
  --out-prefix plots/raw_intake_qc/dfci-g2c.phase1.raw_intake_qc

# Compress & copy initial QC plots to results bucket
tar -czvf data/tarballs/dfci-g2c.phase1.raw_intake_qc.plots.tar.gz plots/raw_intake_qc
gsutil -m cp \
  data/tarballs/dfci-g2c.phase1.raw_intake_qc.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/intake_qc/

# Extract list of QC hard fail samples (missing all phenotype data, reported != inferred sex, etc.)
code/scripts/get_intake_qc_hard_fails.R \
  --qc-tsv data/dfci-g2c.intake_qc.all.tsv.gz \
  --outfile data/dfci-g2c.intake_qc.hard_fail.samples.list

# Build list of global and batch-specific QC exemptions
code/scripts/get_global_qc_exemptions.R \
  --qc-tsv data/dfci-g2c.intake_qc.all.tsv.gz \
  --outfile data/dfci-g2c.intake_qc.global_exemptions.tsv
fgrep -v "wgd_score" \
  data/dfci-g2c.intake_qc.global_exemptions.tsv \
> data/dfci-g2c.intake_qc.batch_exemptions.tsv

# Run batching & QC procedure
zcat data/dfci-g2c.intake_qc.all.tsv.gz > data/dfci-g2c.intake_qc.all.tsv
code/scripts/make_batches.py \
  --match-on batching_sex \
  --match-on batching_pheno \
  --batch-by batching_tissue \
  --batch-by batching_read_length \
  --batch-by wgd_score \
  --batch-by insert_size \
  --batch-by median_coverage \
  --global-qc-cutoffs code/refs/json/dfci-g2c.gatk-sv.global_qc_thresholds.json \
  --global-exemptions data/dfci-g2c.intake_qc.global_exemptions.tsv \
  --batch-qc-cutoffs code/refs/json/dfci-g2c.gatk-sv.batch_qc_thresholds.json \
  --batch-exemptions data/dfci-g2c.intake_qc.batch_exemptions.tsv \
  --custom-qc-fail-samples data/dfci-g2c.intake_qc.hard_fail.samples.list \
  --batch-size 450 \
  --prefix g2c \
  --short-batch-names \
  --outfile data/dfci-g2c.intake_qc.all.post_qc_batching.tsv \
  --batch-names-tsv data/dfci-g2c.gatk-sv.batches.list \
  --batch-membership-tsv data/dfci-g2c.gatk-sv.batch_membership.tsv \
  --fail-reasons-log data/dfci-g2c.intake_qc.all.sample_failure_reasons.tsv \
  --logfile data/dfci-g2c.intake_qc.all.qc_and_batching.log \
  data/dfci-g2c.intake_qc.all.tsv
gzip -f data/dfci-g2c.intake_qc.all.post_qc_batching.tsv

# Copy manifest with updated QC & batching labels to main workspace bucket
gsutil -m cp \
  data/dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/

# Generate QC plots after global QC
code/scripts/plot_intake_qc.R \
  --qc-tsv data/dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
  --pass-column global_qc_pass \
  --out-prefix plots/global_qc_pass/dfci-g2c.phase1.global_qc_pass

# Generate QC plots after batch-specific QC
code/scripts/plot_intake_qc.R \
  --qc-tsv data/dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --out-prefix plots/batch_qc_pass/dfci-g2c.phase1.batch_qc_pass

# Compress & copy global & batch-specific QC plots to results bucket
tar -czvf data/tarballs/dfci-g2c.phase1.global_qc_pass.plots.tar.gz plots/global_qc_pass
tar -czvf data/tarballs/dfci-g2c.phase1.batch_qc_pass.plots.tar.gz plots/batch_qc_pass
gsutil -m cp \
  data/tarballs/dfci-g2c.phase1.global_qc_pass.plots.tar.gz \
  data/tarballs/dfci-g2c.phase1.batch_qc_pass.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/intake_qc/

# Format batch info and stage in central google bucket for GATK-SV
sed '1d' data/dfci-g2c.gatk-sv.batches.list \
> data/gatksv_batch_info/dfci-g2c.gatk-sv.batches.list
while read bid; do
  awk -v bid=$bid '{ if ($1==bid) print $2 }' \
    data/dfci-g2c.gatk-sv.batch_membership.tsv \
  | sort -V | uniq \
  > data/gatksv_batch_info/sample_lists/$bid.samples.list
done < data/gatksv_batch_info/dfci-g2c.gatk-sv.batches.list
code/scripts/evenSplitter.R \
  -S 5 --shuffle \
  data/gatksv_batch_info/dfci-g2c.gatk-sv.batches.list \
  data/gatksv_batch_info/dfci-g2c.gatk-sv.batches.
for i in $( seq 1 5 ); do
  cat data/gatksv_batch_info/dfci-g2c.gatk-sv.batches.$i | sort -V \
  > data/gatksv_batch_info/dfci-g2c.gatk-sv.batches.w$i.list
  rm data/gatksv_batch_info/dfci-g2c.gatk-sv.batches.$i
done
gsutil -m cp -r \
  data/gatksv_batch_info \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/batch_info

# Copy all other batching & QC files to results bucket for reference
gsutil -m cp \
  data/dfci-g2c.intake_qc.hard_fail.samples.list \
  data/dfci-g2c.intake_qc.global_exemptions.tsv \
  data/dfci-g2c.intake_qc.batch_exemptions.tsv \
  data/dfci-g2c.intake_qc.all.sample_failure_reasons.tsv \
  data/dfci-g2c.intake_qc.all.qc_and_batching.log \
  $MAIN_WORKSPACE_BUCKET/results/intake_qc/

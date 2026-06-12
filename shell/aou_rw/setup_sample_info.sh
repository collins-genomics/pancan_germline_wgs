#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell subroutine to stage sample lists and manifests in an AoU workspace

# Note that this code is designed to be run inside the AoU Researcher Workbench
# It also depends on variables set in aou_gatk.sh


# Note2: this code was defined for the AoU RW v1.0, and it may not translate perfectly to AoU RW v2.0 (Verily Pre)
# These will be updated over time as needed for forward compatability, but those
# updates are not guaranteed to preserve reverse compatability with legacy AoU RW v1.0


# Copy sample info to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/data/sample_info/sample_lists ./

# Copy manifest of CRAM/CRAI paths
gsutil -u $GPROJECT -m cp \
  gs://fc-aou-datasets-controlled/v7/wgs/cram/manifest.csv \
  ./cram_manifest.csv

# Write list of cancers ready to be processed
cat << EOF > cancers.list
pancreas
esophagus
stomach
lung
liver
colorectal
kidney
melanoma
prostate
ovary
breast
uterus
cns
bladder
oral
other
male_controls.p1
male_controls.p2
female_controls.p1
female_controls.p2
multiple_sites
pediatric
ufc_cases
ufc_controls
prostate_topup
male_controls.p3
relatives
october_2024_topup_shard1
october_2024_topup_shard2
october_2024_topup_shard3
october_2024_topup_shard4
october_2024_topup_shard5
g2c.phase1.aou.read_metrics.shard1
g2c.phase1.aou.read_metrics.shard2
g2c.phase1.aou.read_metrics.shard3
g2c.phase1.aou.read_metrics.shard4
g2c.phase1.aou.read_metrics.shard5
EOF
# Note: as of 1/8/24, single-sample processing for cancers were divided among 
# multiple workspaces as follows:
# Main workspace: pancreas, esophagus, stomach, lung, liver, kidney, ovary, cns, male_controls.p1, pediatric, october_2024_topup_shard1, g2c.phase1.aou.read_metrics.shard1
# Second workspace: colorectal, bladder, male_controls.p2, multiple_sites, ufc_cases, male_controls.p3, october_2024_topup_shard2, g2c.phase1.aou.read_metrics.shard2
# Third workspace: melanoma, uterus, female_controls.p1, ufc_controls, relatives, october_2024_topup_shard3, g2c.phase1.aou.read_metrics.shard3
# Fourth workspace: prostate, oral, female_controls.p2, october_2024_topup_shard4, g2c.phase1.aou.read_metrics.shard4
# Fifth workspace: breast, other (GCT + testicular), prostate_topup, october_2024_topup_shard5, g2c.phase1.aou.read_metrics.shard5

# Make .tsv mapping person_id, cram path, and crai path for each cancer type
while read cancer; do
  fgrep -wf sample_lists/$cancer.samples.list cram_manifest.csv \
  | sed 's/,/\t/g' | sort -Vk1,1 \
  > data/cram_paths/$cancer.cram_paths.tsv
done < cancers.list

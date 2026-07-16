#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Code to copy all necessary libraries, tools, and other AoU-related code
# to main AoU workspace bucket

# See aou_prep_wdls.sh for WDL, input .json, and other Cromwell files

# Note that this code is designed to be run *locally* (not on Verily)


# Set up local working directory
EXDIR=`pwd`
WRKDIR=`mktemp -d`
cd $WRKDIR

# Clone G2C repo & checkout branch of interest
export g2c_branch=posthoc_filtering
git clone git@github.com:collins-genomics/pancan_germline_wgs.git --branch=$g2c_branch
cd pancan_germline_wgs && \
git rev-parse --short HEAD \
| awk -v FS="\t" '{ print "pancan_germline_wgs", "commit", $1 }' \
> $WRKDIR/libs.version_info.txt && \
cd -

# Clone RLCtools repo & checkout branch of interest
export rlctools_branch=main
git clone git@github.com:RCollins13/RLCtools.git --branch=$rlctools_branch
cd rlctools && \
git rev-parse --short HEAD \
| awk -v FS="\t" '{ print "rlctools", "commit", $1 }' \
>> $WRKDIR/libs.version_info.txt && \
cd -

# Clone GATK-SV repo & checkout release tag of interest
export gatksv_tag=v1.0.1
git clone git@github.com:broadinstitute/gatk-sv.git --branch=$gatksv_tag
echo -e "gatk-sv\ttag\t$gatksv_tag" >> $WRKDIR/libs.version_info.txt

# Make & populate directory of libraries and other tools
for dir in src bin; do
  if [ -e $dir ]; then rm -rf $dir; fi
  mkdir $dir
done
cp -r \
  pancan_germline_wgs/src/g2cpy \
  pancan_germline_wgs/src/G2CR_*.tar.gz \
  RLCtools/RLCtools_*.tar.gz \
  gatk-sv/src/svtk \
  $WRKDIR/src/

# Copy code to main AoU workspace bucket
# Note: must use AoU Google credentials
export rw_bucket=gs://rw-migration-aou-rw-84a0039b
# gcloud auth login
gsutil -m cp -r src $rw_bucket/code/
gsutil -m cp -r \
  pancan_germline_wgs/shell/aou_rw/general_bash_utils.sh \
  pancan_germline_wgs/shell/aou_rw/legacy_rw_1.0 \
  pancan_germline_wgs/shell/aou_rw/aou_bash_utils.sh \
  pancan_germline_wgs/shell/aou_rw/gatksv_bash_utils.sh \
  pancan_germline_wgs/shell/aou_rw/setup_sample_info.sh \
  pancan_germline_wgs/shell/aou_rw/install_packages.sh \
  $rw_bucket/code/refs/
gsutil -m cp \
  pancan_germline_wgs/shell/aou_rw/configure_verily_vm.sh \
  $rw_bucket/code/scripts/
gsutil -m cp \
  pancan_germline_wgs/refs/config/environment.g2c_analysis.yml \
  $rw_bucket/code/refs/config/

# Update commit/tag info for version tracking purposes
gsutil cp $WRKDIR/libs.version_info.txt $rw_bucket/code/refs/


# Clean up
cd $EXDIR
rm -rf $WRKDIR

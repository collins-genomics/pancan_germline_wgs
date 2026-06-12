#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Main variant filtering after callset integration

# Note that this code is designed to run on the All of Us Workbench v2.0 (Verily Pre)
# It may not have reverse compatability with the legacy All of Us Workbench v1.0


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://rw-migration-aou-rw-84a0039b
gcloud storage cp $MAIN_WORKSPACE_BUCKET/code/scripts/configure_verily_vm.sh ./ && \
. configure_verily_vm.sh && \
rm configure_verily_vm.sh


##############################################################
# Compute sample-level global no-call rate per variant class #
##############################################################

# TODO: implement this


################################################
# Curate sample-level covariates for filtering #
################################################

# TODO: implement this


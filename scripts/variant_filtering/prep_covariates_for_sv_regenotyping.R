#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Format sample covariates for SV AD adjustment prior to training SV GT imputation


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Prep covariates for SV GT imputation")
parser$add_argument("--qc-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Intake QC .tsv")
parser$add_argument("--out-tsv", metavar="path", type="character",
                    default="imputation_covariates.tsv",
                    help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.post_qc_batching.tsv.gz",
#              "out_tsv" = "~/scratch/imputation_covariates.test.tsv")

# Load data
df <- load.sample.qc.df(args$qc_tsv)

# Transform variables as needed
df$wgd_plus <- relu(df$wgd_score)
df$wgd_minus <- relu(-df$wgd_score)
df$short_reads <- as.numeric(df$read_length < 150)
df$blood <- as.numeric(df$batching_tissue == "blood")
df$aou <- as.numeric(df$cohort == "aou")
df$g2c_id <- rownames(df)
# TODO: add G2C realignment?

# Retain only covariates needed for SV imputation and write to --out-tsv
df.out <- df[, c("g2c_id", "wgd_plus", "wgd_minus", "insert_size",
                 "mean_coverage", "charr", "chrX_ploidy", "chrY_ploidy",
                 "short_reads", "blood", "aou")]
write.table(df.out, args$out_tsv, col.names=T, row.names=F, quote=F, sep="\t")

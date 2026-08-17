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
G2CR::load.constants("other")


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
parser$add_argument("--hq-samples", metavar="path", type="character",
                    default="hq_samples.list",
                    help="Path to output .txt of high-quality training samples")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz",
#              "out_tsv" = "~/scratch/imputation_covariates.test.tsv",
#              "hq_samples" = "~/scratch/imputation_test.hq_samples.list")

# Load data and subset to samples retained after GATK-HC outlier exclusion
df <- load.sample.qc.df(args$qc_tsv)
df <- df[which(as.logical(df$gatkhc_posthoc_qc_pass)), ]

# Transform variables as needed
df$wgd_plus <- relu(df$wgd_score)
df$wgd_minus <- relu(-df$wgd_score)
df$short_reads <- as.numeric(df$read_length < 150)
df$blood <- as.numeric(df$batching_tissue == "blood")
df$aou <- as.numeric(df$cohort == "aou")
df$g2c_id <- rownames(df)
df$realigned <- as.numeric(as.logical(remap(df$cohort, cohort.realigned, default=FALSE)))

# Retain only covariates needed for SV imputation and write to --out-tsv
df.out <- df[, c("g2c_id", "wgd_plus", "wgd_minus", "insert_size",
                 "mean_coverage", "charr", "chrX_ploidy", "chrY_ploidy",
                 "aou", "realigned")]
write.table(df.out, args$out_tsv, col.names=T, row.names=F, quote=F, sep="\t")

# Define high-quality subset of non-admixed samples for model training
df.hq <- df[which(df$n_grafpop_snps > 58000 & df$n_grafpop_snps < 63000
                  & (df$pct_EUR > 0.8 | df$pct_EUR < 0.3)
                  & (df$pct_AFR > 0.6 | df$pct_AFR < 0.15)
                  & (df$pct_ASN > 0.9 | df$pct_ASN < 0.25)
                  & apply(abs(df[, grep("chr[0-9]+_ploidy", colnames(df))] - 2) < 0.25, 1, all)
                  & df$wgd_score > -0.25 & df$wgd_score < 0.1
                  & df$hq_het_rate > 0.999 & df$mean_ref_ab_hom_alt < 0.002
                  & df$inconsistent_ab_het_rate < 0.05
                  & df$mean_coverage > 20 & df$mean_coverage < 40
                  & df$charr < 0.0075
                  & df$blood == 1
                  & df$short_reads == 0),
            ]


#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper script to filter imputed genotypes by comparing them to existing GTs
# This is much faster than doing serially in pysam


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)



###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Impute GTs for a single SV")
parser$add_argument("--old-gts", metavar=".tsv", type="character", required=TRUE,
                    help="File with variant ID, sample, GT, GQ in current VCF")
parser$add_argument("--imputed-gts", metavar=".tsv", type="character", required=TRUE,
                    help="Output from impute_sv_gts.R")
parser$add_argument("--gq-offset", default=0, type="numeric", metavar="float",
                    help=paste("Difference between imputed GQ and old GQ",
                               "before overwriting old GT"))
parser$add_argument("--out-tsv", metavar="path", type="character",
                    required=TRUE, help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("old_gts" = "~/Downloads/current_gts.tsv.gz",
#              "imputed_gts" = "~/Downloads/dfci-g2c.v1.chr19.final_cleanup_INS_chr19_544.imputation_results.tsv",
#              "gq_offset" = 9,
#              "out_tsv" = "~/scratch/imp_filter.test.tsv")

# Read old & new GTs
old <- read.table(args$old_gts, header=T, check.names=F, sep="\t", comment.char="")
new <- read.table(args$imputed_gts, header=T, check.names=F, sep="\t", comment.char="")
new.cols <- colnames(new)

# Merge old & new GTs
m <- merge(new, old, by=c("#sv_id", "sample"), all=F, sort=F, suffixes=c("", ".old"))
rm(old, new)

# Exclude identical GTs
m <- m[which(m$GT != m$GT.old), ]

# Always overwrite missing GTs
# Otherwise, enforce --gq-offset
m <- m[which(m$GT.old == "./." | (m$GQ > as.numeric(m$GQ.old) + args$gq_offset)), ]

# Write retained updates to --out-tsv
write.table(m[, new.cols], args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)

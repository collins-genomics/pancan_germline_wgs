#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Create sample metadata for G2C genotype & variant filtering


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
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate sample metadata for callset filtering")
parser$add_argument("--qc-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Intake QC .tsv")
parser$add_argument("--samples", metavar=".txt", type="character",
                    help="List of sample IDs to retain [default: keep all]")
parser$add_argument("--out-tsv", metavar="path", type="character",
                    default="G2C.sample_metadata.filtering.tsv",
                    help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.post_qc_batching.tsv.gz",
#              "samples" = NULL,
#              "out_tsv" = "~/scratch/dfci-g2c.filtering_metadata.test.tsv")

# Load data
qc <- load.sample.qc.df(args$qc_tsv)
if(!is.null(args$samples)){
  keep.sids <- intersect(rownames(qc), unique(read.table(args$samples)[, 1]))
  cat(paste("Retained", prettyNum(length(keep.sids), big.mark=","),
            "samples after enforcing --samples\n"))
  qc <- qc[keep.sids, ]
}

# Transform individual features as needed
qc$relu_wgd <- RLCtools::relu(qc$wgd_score)
qc$relu_neg_wgd <- RLCtools::relu(-qc$wgd_score)
qc$short_readlength <- as.integer(qc$read_length <= 150)
qc$blood_wgs <- as.integer(qc$batching_tissue == "blood")
qc$aou <- as.integer(qc$cohort == "aou")

# Extract and standard-normalize all raw features
final.features <- c("mean_coverage", "relu_wgd", "relu_neg_wgd", "short_readlength",
                    "blood_wgs", "insert_size", "aou", "pct_genome_nondiploid",
                    "charr", "mean_ref_ab_hom_alt")
x <- as.data.frame(apply(qc[, final.features], 2, scale, center=T, scale=T))

# Drop singular or undefined features
keep.cidx <- which(!apply(x, 2, function(vals){
  length(table(vals)) == 1 | all(is.na(vals)) | all(is.nan(vals)) | all(is.infinite(vals))
}))
x <- x[, keep.cidx]
x$`#sid` <- rownames(qc)

# Write features to output .tsv
write.table(x[, c("#sid", setdiff(colnames(x), "#sid"))],
            args$out_tsv, col.names=T, row.names=F, quote=F, sep="\t")

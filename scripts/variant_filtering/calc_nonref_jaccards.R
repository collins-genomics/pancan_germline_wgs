#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Compute Jaccard index for variant carriers among a prespecified set of pairs
# This is a subroutine within the larger indel/SV integration workflow


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=FALSE)
require(argparse, quietly=TRUE)


##################
# Data Functions #
##################
# Load genotype matrix
load.gt.matrix <- function(tsv.in){
  # Load genotype matrix and set rownames
  gt <- read.table(tsv.in, header=T, sep="\t", comment.char="", check.names=F)
  rownames(gt) <- gt[, 1]
  gt[, 1] <- NULL
  return(gt)
}

# Compute non-ref sample Jaccard index for a pair of variants
nonref.jaccard <- function(gt, vids){
  df <- data.frame("x1" = as.integer(RLCtools::parse.gt(gt[vids[1], ]) > 0),
                   "x2" = as.integer(RLCtools::parse.gt(gt[vids[2], ]) > 0))
  df <- df[which(complete.cases(df)), ]
  df <- df[which(df$x1 > 0 | df$x2 > 0), ]
  if(nrow(df) == 0){
    return(0)
  }else{
    return(sum(df$x1 == 1 & df$x2 == 1) / nrow(df))
  }
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Compute carrier Jaccards for variant pairs")
parser$add_argument("--genotype-matrix", required=TRUE, metavar=".tsv",
                    help=paste(".tsv of genotypes for all variants.",
                               "First column must be variant ID."))
parser$add_argument("--pairs-tsv", required=TRUE, metavar=".tsv",
                    help=paste(".tsv with at least two columns specifying which",
                               "pairs of variants should be evaluated. Columns",
                               "after the first two will be ignored."))
parser$add_argument("--out-tsv", metavar="path", type="character",
                    help="Destination of output .tsv with Jaccard indexes",
                    default="stdout")
args <- parser$parse_args()

# # DEV:
# args <- list("genotype_matrix" = "~/Downloads/jaccard.input.tsv.gz",
#              "pairs_tsv" = "~/Downloads/candidate_hits.pairs.tsv.gz",
#              "out_tsv" = "~/scratch/jaccard.test.tsv")

# Read variant pairs to process
pairs <- unique(read.table(args$pairs_tsv, header=F, sep="\t")[, 1:2])

# Read genotype matrix
gt <- load.gt.matrix(args$genotype_matrix)

# Compute non-ref Jaccard indexes for all pairs
pairs$jaccard <- sapply(1:nrow(pairs), function(i){
  nonref.jaccard(gt, as.character(pairs[i, ]))
})

# Write Jaccard results to file
colnames(pairs)[1:2] <- c("#vid1", "vid2")
write.table(pairs,
            if(args$out_tsv %in% c("stdout", "-", "/dev/stdout")){""}else{args$out_tsv},
            col.names=T, row.names=F, sep="\t", quote=F)

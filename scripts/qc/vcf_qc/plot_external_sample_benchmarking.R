#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot sample-level genotype benchmarking versus an external dataset


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
load.constants("all")


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot sample-level benchmarking versus external data")
parser$add_argument("--sens-tsv", metavar=".tsv", action="append",
                    help=paste("Sensitivity benchmarking results. Can be",
                               "specified multiple times for multiple",
                               "evaluation interval sets."))
parser$add_argument("--ppv-tsv", metavar=".tsv", action="append",
                    help=paste("PPV benchmarking results. Can be",
                               "specified multiple times for multiple",
                               "evaluation interval sets."))
parser$add_argument("--set-name", metavar="string", action="append",
                    help=paste("Optional names for each set of benchmarking",
                               "results supplied"))
parser$add_argument("--ref-title", metavar="path", type="character",
                    help="String title for comparison dataset/cohort")
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
parser$add_argument("--subset-samples", metavar=".txt",
                    help=paste("Optional flat .txt file with a list of sample IDs",
                               "to keep. By default, all samples will be retained."))
parser$add_argument("--custom-constants", metavar=".R", type="character",
                    help="Optional file of custom constants to use for plotting")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("sens_tsv" = c("~/scratch/sample_gt_bench_plot_dev/giab_easy.external_srwgs.dfci-g2c.v1.initial_qc.chr22.gt_comparison.distrib.merged.tsv.gz",
#                             "~/scratch/sample_gt_bench_plot_dev/giab_hard.external_srwgs.dfci-g2c.v1.initial_qc.chr22.gt_comparison.distrib.merged.tsv.gz"),
#              "ppv_tsv" = c("~/scratch/sample_gt_bench_plot_dev/giab_easy.external_srwgs.dfci-g2c.v1.initial_qc.chr22.gt_comparison.distrib.merged.tsv.gz",
#                            "~/scratch/sample_gt_bench_plot_dev/giab_hard.external_srwgs.dfci-g2c.v1.initial_qc.chr22.gt_comparison.distrib.merged.tsv.gz"),
#              "set_name" = c("Easy", "Hard"),
#              "ref_title" = "external srWGS",
#              "common_af" = 0.001,
#              "subset_samples" = NULL,
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/g2c.qc.test")

# # DEV (single vc):
# args <- list("sens_tsv" = c("~/scratch/dbg_dat/dfci-ufc.sv.v1.initial_qc.external_lrwgs.Easy.ppv_by_freq.merged.tsv.gz",
#                             "~/scratch/dbg_dat/dfci-ufc.sv.v1.initial_qc.external_lrwgs.Hard.ppv_by_freq.merged.tsv.gz"),
#              "ppv_tsv" = c("~/scratch/dbg_dat/dfci-ufc.sv.v1.initial_qc.external_lrwgs.Easy.sensitivity_by_freq.merged.tsv.gz",
#                            "~/scratch/dbg_dat/dfci-ufc.sv.v1.initial_qc.external_lrwgs.Hard.sensitivity_by_freq.merged.tsv.gz"),
#              "set_name" = c("Easy", "Hard"),
#              "ref_title" = "external lrWGS",
#              "common_af" = 0.01,
#              "subset_samples" = NULL,
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/ufc.sv.qc.test")

# Load custom constants if optioned
if(!is.null(args$custom_constants)){
  source(args$custom_constants)
}

# Load list of samples to retain, if optioned
if(!is.null(args$subset_samples)){
  keep.samples <- unique(read.table(args$subset_samples, header=F)[, 1])
}else{
  keep.samples <- NULL
}

# Load sensitivity & PPV data
sens.dat <- load.gt.benchmark.tsvs(args$sens_tsv, args$set_name,
                                   keep.samples=keep.samples)
ppv.dat <- load.gt.benchmark.tsvs(args$ppv_tsv, args$set_name,
                                  keep.samples=keep.samples)

# Get nonredundant list of samples considered in either sensitivity or PPV
samples <- c()
if(!is.null(sens.dat)){
  samples <- unique(c(samples, unlist(lapply(sens.dat, function(d){d$sample}))))
}
if(!is.null(ppv.dat)){
  samples <- unique(c(samples, unlist(lapply(ppv.dat, function(d){d$sample}))))
}
n.samples <- length(samples)

# Define palette for sets
set.colors <- rev(RLCtools::categorical.rainbow(2))

# Prepare collection dataframe for summary stat logging
ss.df <- data.frame("analysis"=character(), "measure"=character(),
                    "value"=numeric(), "n"=numeric())

# Transform reference title to safe string for filenames
out.prefix.base <- paste(args$out_prefix,
                         sub("[ /]+", "_", tolower(args$ref_title)), sep=".")

# Plot sensitivity, if provided
if(length(sens.dat) > 0){
  sens.df <- plot.all.gt.bench.strata(sens.dat,
                                      paste(out.prefix.base, "sensitivity", sep="."),
                                      set.colors,
                                      ref.title=args$ref_title,
                                      metric.name="Sensitivity")
  ss.df <- rbind(ss.df, sens.df)
}

# Plot PPV, if provided
if(length(ppv.dat) > 0){
  ppv.df <- plot.all.gt.bench.strata(ppv.dat, paste(out.prefix.base, "ppv", sep="."),
                                     set.colors,
                                     ref.title=args$ref_title,
                                     metric.name="PPV",
                                     title.preps=c("of", "vs"))
  ss.df <- rbind(ss.df, ppv.df)
}

# Write summary stats to output file for logging
ss.df$n[which(is.na(ss.df$n))] <- n.samples
for(col in colnames(ss.df)){
  ss.df[, col] <- unlist(ss.df[, col])
}
ss.df <- ss.df[with(ss.df, order(analysis, measure)), ]
colnames(ss.df)[1] <- paste("#", colnames(ss.df)[1], sep="")
write.table(ss.df, paste(out.prefix.base, "summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

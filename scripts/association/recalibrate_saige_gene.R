#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins, Noah Fields, and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Recalibrate P-values for each criteria & AF bin from SAIGE-Gene+


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(RLCtools, quietly=TRUE)


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Recalibrate SAIGE-Gene+ statistics")
parser$add_argument("-i", "--input-stats-tsv", required=TRUE, metavar=".tsv",
                    help="Input .tsv of SAIGE-Gene+ statistics")
parser$add_argument("-o", "--output-stats-tsv", required=TRUE, metavar=".tsv",
                    help="Output .tsv of recalibrated SAIGE-Gene+ statistics")
args <- parser$parse_args()

# # DEV
# args <- list("input_stats_tsv" = "~/Downloads/kidney.saige (1).tsv",
#              "output_stats_tsv" = "~/scratch/kidney.saige.recal.tsv")

# Load original statistics
s <- read.table(args$input_stats_tsv, header=T, sep="\t", comment.char="")

# Determine unique combinations of `Group` and `max_MAF`
strata <- unique(s[, c("Group",  "max_MAF")])

# Recalibrate each stratum to account for differences in underlying data density
s.recal <- do.call("rbind", apply(strata, 1, function(sv){
  # Subset statistics to stratum of interest
  g <- as.character(sv[1])
  af <- as.numeric(sv[2])
  sub.df <- s[which(s$Group == g & s$max_MAF == af), ]

  # Backcalculate Z-score from SAIGE-Gene overall P-value and burden beta
  sub.df$raw_Zscore <- qnorm(sub.df$Pvalue/2, lower.tail=(sub.df$BETA_Burden<0))
  spa.df <- RLCtools::saddlepoint.adj(sub.df$raw_Zscore)
  colnames(spa.df) <- c("recalibrated_Zscore", "recalibrated_p")

  # Return dataframe for this stratum extended to include recalibrated statistics
  as.data.frame(cbind(sub.df, spa.df))
}))

# Organize recalibrated results & write to --output-stats-tsv
first.cols <- c("Region", "Group", "max_MAF", "recalibrated_p",
                "BETA_Burden", "SE_Burden")
drop.cols <- colnames(s.recal)[grep("Zscore", colnames(s.recal))]
out.cols <- c(first.cols, setdiff(setdiff(colnames(s.recal), first.cols), drop.cols))
out.df <- s.recal[order(s.recal$recalibrated_p), out.cols]
col.name.map <- c("Region" = "#gene",
                  "Group" = "criteria",
                  "max_MAF" = "max_AF",
                  "BETA_Burden" = "beta",
                  "SE_Burden" = "beta_SE",
                  "Pvalue" = "original_p",
                  "Pvalue_Burden" = "burden_test_p",
                  "Pvalue_SKAT" = "SKAT_p",
                  "Number_rare" = "n_rare_variants",
                  "Number_ultra_rare" = "n_ultra_rare_variants")
colnames(out.df) <- RLCtools::remap(colnames(out.df), col.name.map)
write.table(out.df, args$output_stats_tsv, col.names=T, row.names=F, quote=F, sep="\t")

#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Define outlier samples based on multiples of IQR of variant counts

# Optionally, can further stratify samples based on categorical labels
# (e.g., ancestry) before defining outliers


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
load.constants("colors")
require(RLCtools, quietly=TRUE)


##################
# Data Functions #
##################
# Load variant count data per sample
# Expected input is .tsv format with columns = sample, variant_type, count
load.counts <- function(tsv.in, header=T){
  counts <- read.table(tsv.in, header=header, sep="\t", comment.char="")
  colnames(counts) <- c("sid", "vtype", "k")
  return(counts)
}

# Load optional sample categorical labels
# Expected input is .tsv format with columns = sample, category_label
load.labels <- function(tsv.in, header=T){
  labels <- read.table(tsv.in, header=header, sep="\t", comment.char="")
  l.v <- as.character(labels[, 2])
  names(l.v) <- as.character(labels[, 1])
  return(l.v)
}

# Compute outlier filter thresholds and define outlier samples
define.outliers <- function(vals, n.iqr=6, filter.lower=TRUE){
  # Compute distributional metrics & thresholds
  quarts <- quantile(vals, probs=c(0.25, 0.75))
  iqr <- IQR(vals)
  bw <- n.iqr * iqr
  thresh.upper <- max(quarts) + bw
  thresh.lower <- if(filter.lower){min(quarts) - bw}else{NA}
  out.upper <- names(vals)[which(vals > thresh.upper)]
  out.lower <- if(filter.lower){names(vals)[which(vals < thresh.lower)]}else{c()}
  return(list("iqr.bounds" = quarts,
              "thresholds" = c(thresh.lower, thresh.upper),
              "outliers" = sort(unique(c(out.upper, out.lower)))))
}


######################
# Plotting Functions #
######################
# Visualize distribution of a single count stratum and annotate with outlier cutoffs
plot.count.distrib <- function(vals, thresholds, n.fail, n.samples, pdf.out, title=NULL){
  out.pct.lab <- paste(round(100 * n.fail / n.samples, 2), "%", sep="")
  pdf(pdf.out, height=2, width=2.5)
  density.w.outliers(vals, style="hist", color=DFCI.colors[["darkblue"]],
                     title=title, title.line=0.9, x.title="Variants",
                     x.title.line=-0.1, x.label.units="count",
                     min.x.label.length=2,
                     y.title="Samples", y.title.line=1.5,
                     parmar=c(1.85, 3.4, 1.9, 0.75))
  abline(v=thresholds, lwd=2, col=DFCI.colors[["yellow"]])
  mtext(paste(prettyNum(n.fail, big.mark=","), " sample",
              if(n.fail!=1){"s"}else{""}, " (", out.pct.lab,
              ") failed", sep=""),
        3, cex=5/6, font=3, col=DFCI.colors[["yellow"]], line=0.1)
  dev.off()
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Define outlier samples by variant counts")
parser$add_argument("--counts-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Three-column .tsv of sample id, variant type, count")
parser$add_argument("--sample-labels-tsv", metavar=".tsv", type="character",
                    help=paste("Optional .tsv mapping samples to categories for",
                               "partitioned outlier definition (e.g., ancestry groups)"))
parser$add_argument("--n-iqr", metavar="float", type="numeric", default=6,
                    help="Number of IQR multiples for defining variant count outliers")
parser$add_argument("--no-lower-filter", action="store_true", default=FALSE,
                    help="Disable filtering of outliers below the lower threshold")
parser$add_argument("--plot", action="store_true", default=FALSE,
                    help="Visualize & annotate count distributions")
parser$add_argument("--plot-title-prefix", metavar="string", type="character",
                    help="String prefix for title of all optional plots")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="Path/output prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("counts_tsv" = "~/scratch/sv_count_dbg/dfci-ufc.v1.postCleanupPart1.counts.subsetted.collapsed.tsv",
#              "sample_labels_tsv" = "~/scratch/sv_count_dbg/dfci-g2c.intake_pop_labels.aou_split.tsv",
#              "n_iqr" = 3,
#              "no_lower_filter" = FALSE,
#              "plot" = TRUE,
#              "plot_title_prefix" = "UFC-SV",
#              "out_prefix" = "~/scratch/dfci-ufc.v1.sv")

# Read count data
counts <- load.counts(args$counts_tsv)
n.samples <- length(unique(counts$sid))
vtypes <- sort(unique(counts$vtype))
cat(paste("Loaded", prettyNum(n.samples, big.mark=","), "samples and",
          prettyNum(length(vtypes), big.mark=","), "variant types from",
          "--counts-tsv\n"))

# Annotate with sample categories, if optioned
if(!is.null(args$sample_labels_tsv)){
  labels <- load.labels(args$sample_labels_tsv)
  sctgs <- sort(unique(labels))
  counts <- counts[which(counts$sid %in% names(labels)), ]
  counts$sctg <- labels[counts$sid]
  n.samples <- length(unique(counts$sid))
  cat(paste("Found", prettyNum(n.samples, big.mark=","), "overlapping samples",
            "and", prettyNum(length(sctgs), big.mark=","), "sample categories",
            "in --sample-labels-tsv\n"))
}else{
  sctgs <- NULL
}

# Process each stratum in parallel
all.outliers <- c()
for(vtype in vtypes){
  stratum <- counts[which(counts$vtype == vtype), ]
  stratum.title <- stratum.title.base <- gsub("^ ", "", paste(args$plot_title_prefix, vtype))
  if(is.null(sctgs)){
    s.v <- as.numeric(stratum$k)
    names(s.v) <- stratum$sid
    outlier.data <- define.outliers(s.v, n.iqr=args$n_iqr,
                                    filter.lower=!args$no_lower_filter)
    all.outliers <- sort(unique(c(all.outliers, outlier.data$outliers)))
    if(args$plot){
      n.out <- length(outlier.data$outliers)
      n.samples <- length(unique(stratum$sid))
      pdf.out <- paste(args$out_prefix, vtype, "count_distrib_w_outlier_thresholds",
                       "pdf", sep=".")
      plot.count.distrib(s.v, outlier.data$thresholds, n.out, n.samples,
                         pdf.out=pdf.out, title=stratum.title)
    }
  }else{
    for(sctg in sctgs){
      stratum <- counts[which(counts$vtype == vtype & counts$sctg == sctg), ]
      stratum.title <- paste(stratum.title.base, "in", sctg)
      s.v <- as.numeric(stratum$k)
      names(s.v) <- stratum$sid
      outlier.data <- define.outliers(s.v, n.iqr=args$n_iqr,
                                      filter.lower=!args$no_lower_filter)
      all.outliers <- sort(unique(c(all.outliers, outlier.data$outliers)))
      if(args$plot){
        n.out <- length(outlier.data$outliers)
        n.stratum.samples <- length(unique(stratum$sid))
        pdf.out <- paste(args$out_prefix, vtype, sctg,
                         "count_distrib_w_outlier_thresholds", "pdf", sep=".")
        plot.count.distrib(s.v, outlier.data$thresholds, n.out, n.stratum.samples,
                           pdf.out=pdf.out, title=stratum.title)
      }
    }
  }
}

# Write list of all outliers to text file
n.out.final <- length(all.outliers)
n.out.pct <- round(100 * n.out.final / n.samples, 2)
cat(paste("A total of ", prettyNum(n.out.final, big.mark=","), " samples (",
          n.out.pct, "%) failed any outlier threshold\n", sep=""))
write.table(sort(unique(all.outliers)),
            paste(args$out_prefix, "outliers.samples.list", sep="."),
            col.names=F, row.names=F, quote=F, sep="\t")


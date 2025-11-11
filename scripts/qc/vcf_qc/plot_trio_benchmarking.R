#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot Mendelian genotype benchmarking among complete parent-child trios


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
load.constants("all")


##################
# Plot functions #
##################
# Main wrapper to plot all benchmarking strata
plot.all.trio.bench.strata <- function(bench.dat, out.prefix, set.colors,
                                       denovo.only=FALSE, parmar=c(1, 2.55, 1.8, 2.5)){
  # Check what vcs are present in data
  vcs <- unique(as.character(sapply(bench.dat, function(d){d$class})))

  # Get unique number of samples across all strata
  n.trios <- length(unique(unlist(lapply(bench.dat, function(d){d$family_id}))))

  # Get AF cutoff for common vs. rare
  af.cutoff <- min(unique(clean.breaks(unique(bench.dat[[1]]$freq_bin))))

  # Set mode for analysis from child or full trio POV
  if(denovo.only){
    metric.prefix <- "child_proportion_inherited"
    out.prefix <- paste(out.prefix, "child_pov", sep=".")
    pov <- "child"
    title.prefix <- "Inheritance rate"
    metric.name <- "Proportion inherited"
  }else{
    metric.prefix <- "mendelian_concordance_rate"
    pov <- "trio"
    title.prefix <- "Mendelian assessment"
    metric.name <- "Concordance rate"
  }
  bench.dat <- lapply(bench.dat, function(bdf){
    pass.cidxs <- grep("^pass_", colnames(bdf))
    if(denovo.only){
      bdf$pass <- bdf$pass_transmitted
      other.cols <- setdiff(colnames(bdf)[-c(1:4)], "denovo")
      zero.cols <- other.cols[grep("^pass", other.cols, invert=T)]
      bdf[, zero.cols] <- 0
    }else{
      bdf$pass <- apply(bdf[, pass.cidxs], 1, sum)
    }
    bdf[, pass.cidxs] <- NULL
    return(bdf)
  })

  # Prepare collection dataframe for summary stat logging
  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.prefix <- "trio_concordance"

  # Plot all variants
  plot.dat <- get.gt.bench.plot.data(bench.dat, trio.mode=TRUE)
  ss.df <- rbind(ss.df,
                 calc.gt.bench.ss(plot.dat, paste(ss.prefix, "all", sep="."),
                                  trio.mode=TRUE))
  vpg <- calc.gt.bench.vpg(plot.dat)
  plot.title <- paste(title.prefix, " of\n", clean.numeric.labels(vpg),
                      " variants/", pov, " in ", clean.numeric.labels(n.trios),
                      " trios", sep="")
  pdf(paste(out.prefix, "all.barplot.pdf", sep="."), height=2.25, width=2.2)
  plot.gt.bench(plot.dat=plot.dat,
                strata.names=c("Rare", "Common"),
                set.colors=set.colors,
                metric.name=metric.name,
                title=plot.title, trio.mode=TRUE,
                parmar=parmar)
  dev.off()

  # Plot each variant class individually, if present
  for(vc in vcs){
    plot.dat <- get.gt.bench.plot.data(bench.dat, vc=vc, trio.mode=TRUE)
    ss.df <- rbind(ss.df,
                   calc.gt.bench.ss(plot.dat, paste(ss.prefix, vc, sep="."),
                                    trio.mode=TRUE))
    vpg <- calc.gt.bench.vpg(plot.dat)
    plot.title <- paste(title.prefix, " of\n",
                        clean.numeric.labels(vpg, min.label.length=2, acceptable.decimals=1),
                        " ", gsub("I", "i", paste(var.class.abbrevs[vc], "s", sep="")),
                        " / ", pov, " in ", clean.numeric.labels(n.trios),
                        " trios", sep="")
    pdf(paste(out.prefix, vc, "barplot.pdf", sep="."), height=2.25, width=2.2)
    plot.gt.bench(plot.dat=plot.dat,
                  strata.names=c("Rare", "Common"),
                  set.colors=set.colors,
                  metric.name=metric.name,
                  title=plot.title, trio.mode=TRUE,
                  parmar=parmar)
    dev.off()
  }

  # Plot short variants by subclass, if present
  if(length(intersect(c("snv", "indel"), vcs)) > 1){
    short.vscs <- unique(as.character(sapply(bench.dat, function(d){
      d$subclass[which(d$class %in% c("snv", "indel"))]
    })))

    plot.dat <- lapply(short.vscs, function(vsc){
      calc.gt.bench.stats(bench.dat, vsc=vsc, trio.mode=TRUE)
    })
    names(plot.dat) <- short.vscs
    short.ss.df <- calc.gt.bench.ss(plot.dat, ss.prefix, summarize.all=FALSE, trio.mode=TRUE)
    ss.df <- rbind(ss.df, short.ss.df)
    vpg <- calc.gt.bench.vpg(plot.dat)
    plot.title <- paste(title.prefix, " of\n",
                        clean.numeric.labels(vpg, min.label.length=2, acceptable.decimals=1),
                        " short vars/", pov, " in ", clean.numeric.labels(n.trios),
                        " trios", sep="")
    pdf(paste(out.prefix, "short_variants.by_subclass.barplot.pdf", sep="."),
        height=2.25, width=1.05+(1.3*length(short.vscs)/3))
    plot.gt.bench(plot.dat=plot.dat,
                  strata.names=var.subclass.abbrevs[short.vscs],
                  set.colors=set.colors,
                  metric.name=metric.name,
                  title=plot.title, trio.mode=TRUE,
                  parmar=parmar)
    dev.off()
  }

  # Plot SVs by subclass, if present
  if("sv" %in% vcs){
    sv.vscs <- unique(as.character(sapply(bench.dat, function(d){
      setdiff(d$subclass[which(d$class =="sv")], "CNV")
    })))

    plot.dat <- lapply(sv.vscs, function(vsc){
      calc.gt.bench.stats(bench.dat, vsc=vsc, trio.mode=TRUE)
    })
    names(plot.dat) <- sv.vscs
    sv.ss.df <- calc.gt.bench.ss(plot.dat, ss.prefix, summarize.all=FALSE,
                                 trio.mode=TRUE)
    ss.df <- rbind(ss.df, sv.ss.df)
    vpg <- calc.gt.bench.vpg(plot.dat)
    plot.title <- paste(title.prefix, " of\n",
                        clean.numeric.labels(vpg, min.label.length=2, acceptable.decimals=1),
                        " SVs/", pov, " in ", clean.numeric.labels(n.trios),
                        " trios", sep="")
    pdf(paste(out.prefix, "SVs.by_subclass.barplot.pdf", sep="."),
        height=2.25, width=1.05+(1.4*length(sv.vscs)/3))
    plot.gt.bench(plot.dat=plot.dat,
                  strata.names=var.subclass.abbrevs[sv.vscs],
                  set.colors=set.colors,
                  metric.name=metric.name,
                  title=plot.title, trio.mode=TRUE,
                  parmar=parmar)
    dev.off()
  }

  # Return summary statistics for logging
  med.ridx <- which(ss.df$measure == "median")
  ss.df$measure[med.ridx] <- paste(ss.df$measure[med.ridx], metric.prefix, sep="_")
  ss.df$measure[-med.ridx] <- paste(metric.prefix, ss.df$measure[-med.ridx], sep="_")
  return(ss.df)
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot Mendelian benchmarking from trios")
parser$add_argument("--bench-tsv", metavar=".tsv", action="append",
                    help=paste("Mendelian concordance benchmarking results.",
                               "Can be specified multiple times for multiple",
                               "evaluation interval sets."))
parser$add_argument("--set-name", metavar="string", action="append",
                    help=paste("Optional names for each set of benchmarking",
                               "results supplied"))
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
parser$add_argument("--custom-constants", metavar=".R", type="character",
                    help="Optional file of custom constants to use for plotting")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("bench_tsv" = c("~/Downloads/trio_bench_dev_data/trio_bench.Easy.concordance_distribution.merged.tsv.gz",
#                              "~/Downloads/trio_bench_dev_data/trio_bench.Hard.concordance_distribution.merged.tsv.gz"),
#              "set_name" = c("Easy", "Hard"),
#              "common_af" = 0.001,
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/g2c.qc.test")

# Load custom constants if optioned
if(!is.null(args$custom_constants)){
  source(args$custom_constants)
}

# Load trio concordance data
bench.dat <- load.gt.benchmark.tsvs(args$bench_tsv, args$set_name,
                                    key.cols=1:4, trio.mode=TRUE)

# Get nonredundant list of trio IDs
trios <- unique(unlist(lapply(bench.dat, function(d){d$family_id})))
n.trios <- length(trios)

# Define palette for sets
set.colors <- rev(RLCtools::categorical.rainbow(2))

# Transform reference title to safe string for filenames
out.prefix.base <- paste(args$out_prefix, "trio_concordance", sep=".")

# Plot overall Mendelian concordance
ss.df <- plot.all.trio.bench.strata(bench.dat, out.prefix.base, set.colors)

# Plot proportion of child variants that were inherited from a parent
ss.df.2 <- plot.all.trio.bench.strata(bench.dat, out.prefix.base,
                                      set.colors, denovo.only=TRUE)
ss.df <- rbind(ss.df, ss.df.2)

# Write summary stats to output file for logging
ss.df$n[which(is.na(ss.df$n))] <- n.trios
ss.df <- ss.df[with(ss.df, order(analysis, measure)), ]
colnames(ss.df)[1] <- paste("#", colnames(ss.df)[1], sep="")
write.table(ss.df, paste(out.prefix.base, "summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)


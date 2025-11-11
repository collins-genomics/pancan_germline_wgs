#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot pointwise site-level benchmarking statistics (as opposed to summary
# metrics; for those, see plot_site_benchmarking_summary.R)

# For plotting practicality, this script subsets to common variant sites
# if that is not already done by the user upstream


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
require(viridis, quietly=TRUE)
load.constants("all")


##################
# Data Functions #
##################
# Load a site benchmarking BED file and subset to minimal required columns
read.bed <- function(bed.in, common_af=0){
  numeric.cols <- c("af", "match_af")
  keep.cols <- c("vc", "vsc", numeric.cols)
  df <- read.table(bed.in, header=T, sep="\t", comment.char="",
                   check.names=F, quote="")[, keep.cols]
  df[is.na(df$af), "af"] <- 0
  df[is.na(df$match_af), "match_af"] <- 0
  df[, numeric.cols] <- as.data.frame(apply(df[, numeric.cols], 2, as.numeric))
  df <- df[which(df$af >= common_af | df$match_af >= common_af), ]
  return(as.data.frame(df))
}


######################
# Plotting Functions #
######################
# Custom palette function for AF scatterplot topo maps
af.topo.pal <- function(N){
  c("white", rev(inferno(N-1)))
}

# AF correlation scatterplot
plot.af.scatter <- function(df, title=NULL, common.af=NULL, ref.prefix=NULL,
                            pt.alpha=NULL, pt.cex=NULL, pt.pch=1, legend.wex=0.08,
                            legend.buffer.wex=0.07, e=10e-8, stat.lab.hex=0.04,
                            style="scatter", max.points=1000000,
                            parmar=c(2, 3, 1, 2.4)){
  # Bottom out AFs at floor(log10(common AF)), if optioned
  plot.vals <- df[, c("af", "match_af")]
  plot.vals <- plot.vals[complete.cases(plot.vals), ]

  # Randomly downsample to max.points for scatterplot style
  if(style == "scatter" & nrow(plot.vals) > max.points){
    set.seed(2025)
    plot.vals <- plot.vals[sample(1:nrow(plot.vals), max.points, replace=F), ]
    downsampled <- TRUE
  }else{
    downsampled <- FALSE
  }

  af.d <- apply(plot.vals, 1, function(v){max(v + e) / min(v + e)})
  af.r2 <- cor(plot.vals$af, plot.vals$match_af)^2
  plot.vals <- log10(plot.vals + e)
  min.af <- if(!is.null(common.af)){floor(log10(common.af))}else{min(plot.vals)}
  plot.vals <- apply(plot.vals, 2, function(v){
    v[which(v < min.af)] <- min.af
    return(v)
  })

  # Dynamically determine point properties
  if(style == "scatter"){
    pt.params <- scatterplot.point.params(nrow(plot.vals))
    if(is.null(pt.cex)){
      pt.cex <- pt.params$cex
    }
    if(is.null(pt.pch)){
      pt.pch <- pt.params$pch
    }
    if(is.null(pt.alpha)){
      pt.alpha <- pt.params$alpha
    }
  }else{
    pt.params <- NULL
  }

  # Prep plot area
  legend.x0 <- legend.buffer.wex*abs(min.af)
  legend.x1 <- legend.x0+(legend.wex*abs(min.af))
  prep.plot.area(c(min.af, 0), c(min.af, 0), parmar=parmar, xaxs="r", yaxs="r")

  # Get axis parameters
  ax.tick.at <- 0:min.af
  ax.tick.labs <- paste('"', paste(100 * (10^(0:-3)), "%", sep=""), '"', sep="")
  if(length(ax.tick.at) > length(ax.tick.labs)){
    ax.tick.labs <- c(ax.tick.labs,
                      paste("10 ^ -", length(ax.tick.labs):length(ax.tick.at)), sep="")
  }else{
    ax.tick.labs <- head(ax.tick.labs, length(ax.tick.at))
  }
  ax.tick.labs[length(ax.tick.labs)] <- paste("'' <", ax.tick.labs[length(ax.tick.labs)])

  # Add X axis
  clean.axis(1, at=ax.tick.at, labels=ax.tick.labs, parse.labels=TRUE,
             label.line=-0.85, title="AF (this cohort)", title.line=0)

  # Add Y axis
  if(is.null(ref.prefix)){
    y.title <- "AF (Ref.)"
  }else{
    y.title <- paste("AF (", ref.prefix, ")", sep="")
  }
  clean.axis(2, at=ax.tick.at, labels=ax.tick.labs, parse.labels=TRUE,
             label.line=-0.75, title=y.title, title.line=1.05)

  # Add points
  if(is.null(pt.params)){
    density.topomap(plot.vals[, 1], jitter(plot.vals[, 2]),
                    xlims=c(min.af, 0), ylims=c(min.af, 0),
                    palette.fxn=af.topo.pal,
                    contour.levels=6)
  }else{
    pw.pal <- inferno(101, alpha=pt.alpha, begin=0.075, end=0.9)
    pw.col.idx <- ceiling(100 * log10(sapply(af.d, function(v){min(c(v, 10))}))) + 1
    points(plot.vals, cex=pt.cex, pch=pt.pch, col=pw.pal[pw.col.idx])
  }

  # Add title & annotations
  abline(0, 1, lty=3, col=adjustcolor("white", alpha=0.4))
  mtext(3, text=title, line=0, xpd=T)
  n.vars <- clean.numeric.labels(nrow(plot.vals))
  n.lab <- paste("N=", n.vars, sep="")
  if(downsampled){
    n.lab <- paste(n.lab, "(downsampled)", sep="\n")
  }
  text(x=0, offset=0.2,
       y=par("usr")[3]+((if(downsampled){3}else{2})*stat.lab.hex*diff(par("usr")[3:4])),
       pos=2, cex=5/6, labels=n.lab, xpd=T)
  text(x=min.af, offset=0.2,
       y=par("usr")[4]-(2.5*stat.lab.hex*diff(par("usr")[3:4])),
       pos=4, cex=5/6,
       labels=bquote(italic(R)^2*"="*.(formatC(round(af.r2, 2), digits=2))))

  # Add legend
  if(!is.null(pt.params)){
    legend.y <- quantile(par("usr")[1:2], probs=seq(0.15, 0.8, length.out=101))
    text(x=legend.x1, y=quantile(legend.y, probs=0.965), pos=3, xpd=T,
         labels=bquote(Phi[AF]))
    rect(xleft=legend.x0, xright=legend.x1,
         ybottom=legend.y[-101], ytop=legend.y[-1],
         border=NA, bty="n", col=adjustcolor(pw.pal[-1], alpha.f=10e10), xpd=T)
    rect(xleft=legend.x0, xright=legend.x1,
         ybottom=min(legend.y), ytop=max(legend.y),
         col=NA, xpd=T)
    text(x=mean(c(legend.x0, legend.x1)),
         y=quantile(legend.y, probs=c(0.025, 0.2, 0.4, 0.6, 0.8, 0.95)),
         cex=4.5/6, pos=4, labels=c("~", "2x", "4x", "6x", "8x", ">10x"), xpd=T)
  }

  return(c(af.r2, nrow(plot.vals)))
}


# Plot wrapper function for all site-level pointwise plots
pointwise.plots <- function(df, out.prefix, fname.suffix="all",
                            title="All variants", common.af=0.01,
                            ref.prefix=NULL){

  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.prefix <- gsub("[ ]+", "_", paste(tolower(ref.prefix), fname.suffix, sep="."))

  # AF correlation scatterplot plot as .png
  png(paste(out.prefix, fname.suffix, "af_cor.png", sep="."),
      height=2.25*300, width=2.6*300, res=300)
  m.tmp <- plot.af.scatter(df, title=title, common.af=common.af,
                           ref.prefix=ref.prefix)
  dev.off()

  # Update sumstats
  ss.df[1, ] <- c(paste(ss.prefix, "common_af_cor", sep="."), "r2", m.tmp)

  # AF correlation plot as topo .pdf
  pdf(paste(out.prefix, fname.suffix, "af_cor.topo.pdf", sep="."),
      height=2.25, width=2.3)
  m.tmp <- plot.af.scatter(df, title=title, common.af=common.af,
                           ref.prefix=ref.prefix, style="topo",
                           parmar=c(2, 3, 1, 0.75))
  dev.off()

  return(ss.df)
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot pointwise variant metrics")
parser$add_argument("--snvs", metavar=".bed", type="character",
                    help="SNV site benchmarking BED from BenchmarkSites.wdl")
parser$add_argument("--indels", metavar=".bed", type="character",
                    help="Indel site benchmarking BED from BenchmarkSites.wdl")
parser$add_argument("--svs", metavar=".bed", type="character",
                    help="SV site benchmarking BED from BenchmarkSites.wdl")
parser$add_argument("--combine", action="store_true", default=FALSE,
                    help="Also generate a combined set of plots for all variant types")
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
parser$add_argument("--ref-title", metavar="string", type="character",
                    help="Name of reference cohort (optional)")
parser$add_argument("--custom-constants", metavar=".R", type="character",
                    help="Optional file of custom constants to use for plotting")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
parser$add_argument("--set-name", metavar="string",
                    help=paste("Optional suffix for summary statistics reporting"))
args <- parser$parse_args()

# # DEV:
# args <- list("snvs" = "~/scratch/site_benchmarking_dev/giab_easy.dfci-g2c.v1.gatkhc.initial_qc.chr19_vs_gnomad_v4.common_sites.snvs.bed.gz",
#              "indels" = "~/scratch/site_benchmarking_dev/giab_easy.dfci-g2c.v1.gatkhc.initial_qc.chr19_vs_gnomad_v4.common_sites.indels.bed.gz",
#              "svs" = "~/scratch/site_benchmarking_dev/giab_easy.dfci-g2c.v1.gatkhc.initial_qc.chr19_vs_gnomad_v4.common_sites.svs.bed.gz",
#              "combine" = TRUE,
#              "common_af" = 0.001,
#              "ref_title" = "gnomAD v4.1",
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/qc.test",
#              "set_name" = "easy")

# Ensure at least one of --snvs, --indels, or --svs is present
if(is.null(args$snvs) & is.null(args$indels) & is.null(args$svs)){
  stop("At least one of --snvs, --indels, or --svs must be provided")
}

# Load custom constants if optioned
if(!is.null(args$custom_constants)){
  source(args$custom_constants)
}

# Load & plot SNVs, if provided
if(!is.null(args$snvs)){
  # Load SNV data
  snv.df <- read.bed(args$snvs, args$common_af)

  # Plot SNV metrics
  snv.ss <- pointwise.plots(snv.df, args$out_prefix, fname.suffix="snv",
                            title="Common SNVs", common.af=args$common_af,
                            ref.prefix=args$ref_title)
}else{
  snv.df <- NULL
  snv.ss <- NULL
}

# Load & plot indels, if provided
if(!is.null(args$indels)){
  # Load indel data
  indel.df <- read.bed(args$indels, args$common_af)

  # Plot indel metrics
  indel.ss <- pointwise.plots(indel.df, args$out_prefix, fname.suffix="indel",
                              title="Common indels", common.af=args$common_af,
                              ref.prefix=args$ref_title)

}else{
  indel.df <- NULL
  indel.ss <- NULL
}

# Load & plot SVs, if provided
if(!is.null(args$svs)){
  # Load SV data
  sv.df <- read.bed(args$svs, args$common_af)

  # Plot SV metrics
  sv.ss <- pointwise.plots(sv.df, args$out_prefix, fname.suffix="sv",
                           title="Common SVs", common.af=args$common_af,
                           ref.prefix=args$ref_title)
}else{
  sv.df <- NULL
  sv.ss <- NULL
}

# Combine & plot all variant types, if optioned
if(args$combine){
  # Merge all data
  all.df <- do.call("rbind", list(snv.df, indel.df, sv.df))

  # Plot all metrics
  all.ss <- pointwise.plots(all.df, args$out_prefix, fname.suffix="all",
                            title="All common variants", common.af=args$common_af,
                            ref.prefix=args$ref_title)
}else{
  all.ss <- NULL
}

# Combine summary statistics and write to .tsv
ss.out <- do.call("rbind", list(snv.ss, indel.ss, sv.ss, all.ss))
if(!is.null(args$set_name)){
  if(args$set_name != "All"){
    ss.out$analysis <- paste(ss.out$analysis, tolower(args$set_name), sep=".")
  }
}
colnames(ss.out)[1] <- paste("#", colnames(ss.out)[1], sep="")
write.table(ss.out, paste(args$out_prefix, "summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)


#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot pointwise variant site-level statistics (as opposed to summary metrics;
# for those, see plot_site_summary_metrics.R)

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
# Load a site metrics BED file and subset to minimal required columns
read.bed <- function(bed.in, common_af=0, autosomes.only=T){
  df <- read.table(bed.in, header=T, sep="\t", comment.char="",
                   check.names=F, quote="")
  if(autosomes.only){
    df <- df[df[, 1] %in% c(1:22, paste("chr", c(1:22), sep="")), ]
  }
  keep.cols <- c("af", "freq_het", "freq_hom", "hwe")
  df[, keep.cols] <- as.data.frame(apply(df[, keep.cols], 2, as.numeric))
  df <- df[which(df$af >= common_af), c("vid", keep.cols)]
  df <- df[complete.cases(df), ]
  df[, c("hwe.x", "hwe.y")] <- t(apply(df[, c("freq_het", "freq_hom")], 1,
                                       function(freqs){calc.hwe.xy(freqs[1], freqs[2])}))
  rownames(df) <- df$vid; df$vid <- NULL
  return(as.data.frame(df))
}


######################
# Plotting Functions #
######################
# Custom palette function for HWE topo maps
hwe.topo.pal <- function(N){
  c("black", viridis(N-1))
}

# HWE ternary plot
hwe.plot <- function(df, title="All variants", pt.cex=NULL, pt.pch=NULL,
                     pt.alpha=NULL, style="scatter", max.points=1000000,
                     parmar=c(1.75, 2, 1.25, 1)){
  # Clean input data and determine overall statistics
  df <- df[which(!is.na(df$hwe.x) & !is.na(df$hwe.y) & !is.na(df$hwe)), ]
  bonf.p <- 0.05 / nrow(df)
  n.all <- nrow(df)
  n.pass <- sum(df$hwe >= bonf.p)
  pct.pass <- n.pass / n.all

  # Downsample scatterplot to max.points if more complete cases are provided
  if(style == "scatter" & nrow(df) > max.points){
    set.seed(2025)
    df <- df[sample(1:nrow(df), max.points), ]
  }

  # Get plot parameters
  ytop <- sin(60 * pi / 180)
  n.p.bins <- floor(10 * -log10(bonf.p))

  # Dynamically determine point properties
  if(style == "scatter"){
    pt.params <- scatterplot.point.params(nrow(df), cex.start=0.3)
    if(is.null(pt.cex)){
      pt.cex <- pt.params$cex
    }
    if(is.null(pt.pch)){
      pt.pch <- pt.params$pch
    }
    if(is.null(pt.alpha)){
      pt.alpha <- pt.params$alpha
    }
    p.pal <- rev(viridis(n.p.bins + 1, begin=0.1, end=0.95))
    pt.cols <- p.pal[sapply(floor(10 * -log10(df$hwe)),
                            function(v){min(v, n.p.bins)}) + 1]
    pt.cols <- adjustcolor(pt.cols, alpha=pt.alpha)
  }else{
    pt.params <- NULL
  }

  # Prep plot area
  prep.plot.area(c(0, 1), c(0, 1), parmar)

  # Add heatmap for topo style
  if(style != "scatter"){
    density.topomap(df$hwe.x, df$hwe.y, contour.levels=6, res=200,
                    palette.fxn=hwe.topo.pal,
                    xlims=c(0, 1), ylims=c(0, ytop))
    polygon(x=c(0, 0.5, 0, 0), y=c(0, ytop, ytop, 0),
            col="white", border="white")
    polygon(x=c(0.5, 1, 1, 0.5), y=c(ytop, ytop, 0, ytop),
            col="white", border="white")
  }
  segments(x0=c(0, 0.5, 1), x1=c(0.5, 1, 0),
           y0=c(0, ytop, 0), y1=c(ytop, 0, 0),
           col=annotation.color, xpd=T)

  # Add bottom X axis
  clean.axis(1, at=seq(0, 1, 0.25), label.units="percent", label.line=-1,
             title="Homozygotes", title.line=-0.25)

  # Add left diagonal X axis
  segments(x0=0, x1=0.5, y0=0, y1=ytop)
  diag.ax.spacing  <- seq(0, 1, 0.25)
  diag.ax.x.at <- diag.ax.spacing * cos(60 * pi / 180)
  diag.ax.y.at <- diag.ax.spacing * sin(60 * pi / 180)
  diag.ax.tck <- 0.025
  diag.ax.x.adj <- diag.ax.tck * cos(30 * pi / 180)
  diag.ax.y.adj <- diag.ax.tck * sin(30 * pi / 180)
  segments(x0=diag.ax.x.at, x1=diag.ax.x.at - diag.ax.x.adj,
           y0=diag.ax.y.at, y1=diag.ax.y.at + diag.ax.y.adj, xpd=T)
  text(x=diag.ax.x.at+0.025, y=diag.ax.y.at+0.015, pos=2, cex=4.5/6,
       labels=paste(seq(0, 100, 25), "%", sep=""), xpd=T)
  text(x=mean(diag.ax.x.at) - (10*diag.ax.x.adj),
       y=mean(diag.ax.y.at) + (10*diag.ax.y.adj),
       srt=60, labels="Heterozygotes", xpd=T)

  # Add legend
  if(style == "scatter"){
    legend.left <- 0.75
    legend.right <- 0.8
    legend.y.breaks <- seq(2*ytop/3, ytop, length.out=4)
    legend.col.breaks <- c(0, floor(10 * -log10(0.05)),
                           floor(10 * -log10(bonf.p)), length(p.pal))
    hwe.leg.labs <- list(format.pval(0.05, equality=">="),
                         format.pval(0.05, equality="<"),
                         format.pval(bonf.p, equality="<"))
    text(x=legend.left+0.025, y=mean(legend.y.breaks[3:4]),
         pos=2, labels=bquote(italic(P)))
    sapply(1:3, function(i){
      k.range <- legend.col.breaks[i]:legend.col.breaks[i+1]
      k.rect.y.breaks <- seq(legend.y.breaks[i], legend.y.breaks[i+1],
                             length.out=length(k.range))
      rect(xleft=rep(legend.left, length(k.range)),
           xright=rep(legend.right, length(k.range)),
           ybottom=k.rect.y.breaks[-length(k.rect.y.breaks)],
           ytop=k.rect.y.breaks[-c(1)], col=p.pal[k.range],
           border=p.pal[k.range], lwd=0.5)
      text(x=legend.left-0.015,
           y=mean(legend.y.breaks[c(i, i+1)]),
           pos=4, cex=5/6, xpd=T,
           labels=if(i==1){
             expression("" > 0.05)
           }else if(i==2){
             expression("" <= 0.05)
           }else if(i==3){
             expression("" < "Bonf.")
           })
    })
    rect(xleft=rep(legend.left, 3), xright=rep(legend.right, 3),
         ybottom=legend.y.breaks[1:3], ytop=legend.y.breaks[2:4],
         col=NA, xpd=T, border="white")
    rect(xleft=legend.left, xright=legend.right,
         ybottom=legend.y.breaks[1], ytop=legend.y.breaks[4],
         col=NA)
  }

  # Add title & subtitle
  mtext(3, line=0.4, text=title)
  n.formatted <- clean.numeric.labels(c(n.pass, n.all), acceptable.decimals=2)
  n.formatted[1] <- gsub("k|M|B|T", "", n.formatted[1])
  subtitle <- paste(n.formatted[1], " / ",
                    n.formatted[2], " (", round(100*pct.pass, 1),
                    "%) pass HWE", sep="")
  mtext(3, cex=5/6, text=subtitle, line=-0.4)

  # Add points
  if(style == "scatter"){
    points(df$hwe.x, df$hwe.y, pch=1, cex=pt.cex, col=pt.cols, xpd=T)
  }

  return(c(pct.pass, n.all))
}

# Generate a plot of peak LD R2 values per variant
ld.plot <- function(df, ld, vc2, title, ld.cutoffs=c(0.2, 0.5, 0.8),
                    ld.bin.labels=c("None", "Weak", "Med.", "Strong"),
                    margin.n.bins=20, right.margin.hex=0.2,
                    bar.buffer.hex=0.015, margin.tick.hex=0.05,
                    parmar=c(2, 2.35, 1, 1)){
  # Retain complete data only
  ld <- ld[which(ld$vid %in% rownames(df)), ]
  if(vc2 == "any"){
    ld <- ld[order(-ld$ld_r2), ]
    ld <- ld[!duplicated(ld$vid), ]
  }else{
    ld <- ld[which(ld$other_vc == vc2), ]
  }
  rownames(ld) <- ld$vid
  plot.df <- merge(df, ld, all=F, sort=F, by=0)[, c("af", "ld_r2")]
  plot.df$af <- log10(plot.df$af)
  plot.df <- plot.df[complete.cases(plot.df), ]

  # Get plot parameters
  min.af <- round(min(plot.df$af), 0)
  bar.pal <- greyscale.palette(margin.n.bins, buffer=margin.n.bins/2)
  right.bar.xmax <- abs(min.af) * right.margin.hex
  bar.x.add <- abs(min.af) * bar.buffer.hex
  margin.tick.len <- abs(min.af) * margin.tick.hex
  ld.bin.breaks <- c(0, ld.cutoffs, 1)

  # Prep plot area
  prep.plot.area(c(min.af, right.bar.xmax+bar.x.add), c(0, 1),
                 parmar=parmar, xaxs="r", yaxs="r")
  segments(x0=par("usr")[1], x1=0, y0=ld.cutoffs, y1=ld.cutoffs, col=annotation.color)

  # Plot points
  density.topomap(plot.df$af, plot.df$ld_r2, xlims=c(min.af, 0), ylims=c(0, 1), ncol=48)

  # Add X axis for AF
  ax.tick.at <- 0:min.af
  ax.tick.labs <- paste('"', paste(100 * (10^(0:-3)), "%", sep=""), '"', sep="")
  if(length(ax.tick.at) > length(ax.tick.labs)){
    ax.tick.labs <- c(ax.tick.labs,
                      paste("10 ^ -",
                            length(ax.tick.labs):length(ax.tick.at)), sep="")
  }else{
    ax.tick.labs <- head(ax.tick.labs, length(ax.tick.at))
  }
  ax.tick.labs[length(ax.tick.labs)] <-
    paste("'' <", ax.tick.labs[length(ax.tick.labs)])
  clean.axis(1, at=ax.tick.at, labels=ax.tick.labs, parse.labels=TRUE,
             label.line=-0.85)
  axis(1, at=mean(c(min.af, 0)), tick=F, line=0, labels=paste(title, "AF"))

  # Add axes & title
  segments(x0=par("usr")[1], x1=0, y0=ld.cutoffs, y1=ld.cutoffs,
           col=adjustcolor("white", alpha=0.3), lty=3)
  if(vc2 == "any"){
    clean.axis(2, at=ld.bin.breaks, label.line=-0.75, title.line=0.1,
               title=bquote("Best" ~ R^2))
    axis(3, at=mean(c(min.af, 0)), labels=paste(title, "LD"), line=-1, tick=F)
  }else{
    vc2.lab <- gsub("I", "i", var.class.abbrevs[vc2])
    clean.axis(2, at=ld.bin.breaks, label.line=-0.75, title.line=0.1,
               title=bquote("Best" ~ R^2 ~ "vs. any" ~ .(vc2.lab)))
    axis(3, at=mean(c(min.af, 0)), tick=F, line=-1,
         labels=paste(title, " LD vs. ", vc2.lab, "s", sep=""))
  }

  # Add marginal histogram on right Y axis
  h <- hist(floor(margin.n.bins * plot.df$ld_r2),
            plot=F, breaks=0:margin.n.bins)$counts
  h.x <- h * right.bar.xmax / max(h)
  rect(xleft=bar.x.add, xright=h.x+bar.x.add,
       ybottom=((1:margin.n.bins)-1)/margin.n.bins,
       ytop=(1:margin.n.bins)/margin.n.bins,
       border=bar.pal, col=bar.pal)
  segments(x0=bar.x.add, x1=bar.x.add+margin.tick.len,
           y0=ld.bin.breaks, y1=ld.bin.breaks)
  segments(x0=bar.x.add, x1=bar.x.add, y0=0, y1=1, lwd=2/3)
  pt.bin.mem <- sapply(plot.df$ld_r2,
                       function(v){which.max(v <= ld.bin.breaks)-1})
  sapply(1:length(ld.bin.breaks), function(i){
    bin.pct <- sum(pt.bin.mem == i) / nrow(plot.df)
    text(x=bar.x.add, y=mean(ld.bin.breaks[i+(0:1)])+(0.03*diff(par("usr")[3:4])),
         cex=5/6, labels=ld.bin.labels[i], xpd=T, pos=4)
    text(x=bar.x.add, y=mean(ld.bin.breaks[i+(0:1)])-(0.045*diff(par("usr")[3:4])),
         cex=4/6, xpd=T, pos=4,
         labels=paste("(", round(100 * bin.pct, 1), "%)", sep=""))
  })

  return(c(mean(plot.df$ld_r2), nrow(plot.df)))
}

# Plot wrapper function for all site-level pointwise plots
pointwise.plots <- function(df, ld, out.prefix, fname.suffix="all",
                            title="All variant"){

  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.prefix <- gsub("[ ]+", "_", fname.suffix)

  # HWE scatterplot as .png
  png(paste(out.prefix, fname.suffix, "hwe.png", sep="."),
      height=2.25*300, width=2.25*300, res=300)
  m.tmp <- hwe.plot(df, title=paste(title, "s", sep=""))
  dev.off()
  ss.df[1, ] <- c(paste(ss.prefix, "common_hwe", sep="."), "pct_pass", m.tmp)

  # HWE topo heatmap as .pdf
  pdf(paste(out.prefix, fname.suffix, "hwe.topo.pdf", sep="."), height=2.25, width=2.25)
  m.tmp <- hwe.plot(df, title=paste(title, "s", sep=""), style="topo")
  dev.off()

  # Exit now if LD stats are not provided
  if(is.null(ld)){
    return(ss.df)
  }

  # Compute LD tag rate
  ld.sub <- ld[which(ld$vid %in% rownames(df)), ]
  n.tagged <- length(unique(ld.sub[which(ld.sub$ld_r2 >= 0.5), "vid"]))
  n.elig <- length(unique(ld.sub$vid))
  ss.df[nrow(ss.df)+1, ] <- c(paste(ss.prefix, "common_ld", "any", sep="."),
                              "tag_rate", n.tagged / n.elig, n.elig)

  # Peak LD vs. AF as topo .pdf (one for each other variant class)
  for(vc2 in c("any", unique(ld$other_vc))){
    pdf(paste(out.prefix, fname.suffix, vc2, "peak_ld.pdf", sep="."),
        height=2.25, width=2.6)
    m.tmp <- ld.plot(df, ld, vc2, title=title)
    dev.off()
    ss.df[nrow(ss.df)+1, ] <- c(paste(ss.prefix, "common_ld", vc2, sep="."),
                                "mean_peak_r2", m.tmp)
  }

  return(ss.df)
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot pointwise variant metrics")
parser$add_argument("--snvs", metavar=".bed", type="character",
                    help="SNV site metrics from clean_site_metrics.py")
parser$add_argument("--indels", metavar=".bed", type="character",
                    help="Indel site metrics from clean_site_metrics.py")
parser$add_argument("--svs", metavar=".bed", type="character",
                    help="SV site metrics from clean_site_metrics.py")
parser$add_argument("--combine", action="store_true", default=FALSE,
                    help="Also generate a combined set of plots for all variant types")
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
parser$add_argument("--ld-stats", metavar=".tsv", type="character",
                    help="Peak LD R2 per variant for each other variant class")
parser$add_argument("--custom-constants", metavar=".R", type="character",
                    help="Optional file of custom constants to use for plotting")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("snvs" = "~/scratch/dfci-g2c.v1.chr22.0.norm.posthoc_filtered.sites.snv.sites.common.bed.gz",
#              "indels" = "~/scratch/dfci-g2c.v1.chr22.0.norm.posthoc_filtered.sites.indel.sites.common.bed.gz",
#              "svs" = "~/scratch/dfci-g2c.v1.chr22.0.norm.posthoc_filtered.sites.sv.sites.common.bed.gz",
#              "combine" = TRUE,
#              "common_af" = 0.001,
#              "ld_stats" = "~/scratch/renamed.common.peak_ld_by_vc.tsv.gz",
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/qc.test")

# Ensure at least one of --snvs, --indels, or --svs is present
if(is.null(args$snvs) & is.null(args$indels) & is.null(args$svs)){
  stop("At least one of --snvs, --indels, or --svs must be provided")
}

# Load custom constants if optioned
if(!is.null(args$custom_constants)){
  source(args$custom_constants)
}

# Load LD stats, if provided
if(!is.null(args$ld_stats)){
  ld <- read.table(args$ld_stats, header=T, sep="\t", check.names=F,
                   quote="", comment.char="")
  colnames(ld)[1] <- gsub("#", "", colnames(ld)[1])
}else{
  ld <- NULL
}

# Load & plot SNVs, if provided
if(!is.null(args$snvs)){
  # Load SNV data
  snv.df <- read.bed(args$snvs, args$common_af)

  # Plot SNV metrics
  snv.ss <- pointwise.plots(snv.df, ld, args$out_prefix, fname.suffix="snv",
                            title="Common SNV")
}else{
  snv.df <- NULL
  snv.ss <- NULL
}

# Load & plot indels, if provided
if(!is.null(args$indels)){
  # Load indel data
  indel.df <- read.bed(args$indels, args$common_af)

  # Plot indel metrics
  indel.ss <- pointwise.plots(indel.df, ld, args$out_prefix, fname.suffix="indel",
                              title="Common indel")

}else{
  indel.df <- NULL
  indel.ss <- NULL
}

# Load & plot SVs, if provided
if(!is.null(args$svs)){
  # Load SV data
  sv.df <- read.bed(args$svs, args$common_af)

  # Plot SV metrics
  sv.ss <- pointwise.plots(sv.df, ld, args$out_prefix, fname.suffix="sv",
                           title="Common SV")
}else{
  sv.df <- NULL
  sv.ss <- NULL
}

# Combine & plot all variant types, if optioned
if(args$combine){
  # Merge all data
  all.df <- do.call("rbind", list(snv.df, indel.df, sv.df))

  # Plot all metrics
  all.ss <- pointwise.plots(all.df, ld, args$out_prefix, fname.suffix="all",
                            title="All common variant")
}else{
  all.ss <- NULL
}

# Combine summary statistics and write to .tsv
ss.out <- do.call("rbind", list(snv.ss, indel.ss, sv.ss, all.ss))
colnames(ss.out)[1] <- paste("#", colnames(ss.out)[1], sep="")
write.table(ss.out, paste(args$out_prefix, "pointwise_summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)


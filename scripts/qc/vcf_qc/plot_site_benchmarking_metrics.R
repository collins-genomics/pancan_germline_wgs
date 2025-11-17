#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot site benchmarking summary metrics


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
# Load & postprocess one or more benchmarking summary .tsvs
load.benchmark.tsvs <- function(tsvs.in, names=NULL){
  res <- lapply(tsvs.in, function(tsv.in){
    bd <- read.compressed.distrib(tsv.in, key.cols=1:4)
    bd$df[, -(1:3)] <- apply(bd$df[, -(1:3)], 2, as.numeric)
    bd$n <- bd$df[, 1:3]
    bd$n$n <- apply(bd$df[, -(1:3)], 1, sum, na.rm=T)
    return(bd)
  })
  names(res) <- if(is.null(names)){1:length(res)}else{names}
  return(res)
}

# Simplify results for a single benchmarking result
extract.bench.res <- function(bd, vc.elig=names(var.class.abbrevs),
                              vsc.elig=names(var.subclass.abbrevs),
                              common.af=0.01, summarize.by.vsc=FALSE){
  # Subset and compress count dataframe
  drop.cols <- if(summarize.by.vsc){c("class", "af")}else{c("class", "subclass")}
  if(summarize.by.vsc){
    b.df <- bd$df[which(bd$df$class %in% vc.elig
                        & clean.breaks(bd$df$af) >= common.af),
                  setdiff(colnames(bd$df), c("class", "af"))]
  }else{
    b.df <- bd$df[which(bd$df$class %in% vc.elig
                        & bd$df$subclass %in% vsc.elig),
                  setdiff(colnames(bd$df), drop.cols)]
  }
  strata <- unique(b.df[, 1])
  norm.df <- as.data.frame(do.call("rbind", lapply(strata, function(sx){
    s.v <- apply(b.df[which(b.df[, 1] == sx), -(1)], 2, sum)
    s.v.n <- s.v / sum(s.v)
    s.v.n[setdiff(names(s.v.n), "no_match")]
  })))
  rownames(norm.df) <- strata

  # Get corresponding sample sizes
  if(summarize.by.vsc){
    b.n <- bd$n[which(bd$n$class %in% vc.elig
                      & clean.breaks(bd$n$af) >= common.af),
                setdiff(colnames(bd$n), drop.cols)]
  }else{
    b.n <- bd$n[which(bd$n$class %in% vc.elig
                      & bd$n$subclass %in% vsc.elig),
                setdiff(colnames(bd$n), drop.cols)]
  }
  b.n <- sapply(strata, function(sx){sum(b.n[which(b.n[, 1] == sx), "n"])})
  return(list("df" = norm.df, "n" = b.n, breaks=bd$breaks))
}

# Compress benchmarking data into a single summary value
compress.benchmarking <- function(bd, vc.elig=names(var.class.abbrevs),
                                  vsc.elig=names(var.subclass.abbrevs),
                                  min.af=0){
  # Filter input data
  df <- bd$df[which(bd$df$class %in% vc.elig
                    & bd$df$subclass %in% vsc.elig), ]
  if("af" %in% colnames(df)){
    df <- df[which(clean.breaks(df$af) > min.af), ]
  }
  df <- df[, -(1:3)]

  # Count total variants and those with no match
  n.all <- sum(df, na.rm=T)
  n.miss <- sum(df[, 1], na.rm=T)
  n.hit <- n.all - n.miss

  # Return hit rate & total number of variants
  if(n.all > 0){
    c(n.hit / n.all, n.all)
  }else{
    c(NA, n.all)
  }
}

# Collect all common benchmarking summary statistics for a single type of inputs
gather.bench.sumstats <- function(bench.dat, ref.title="ref. cohort",
                                  metric="concordance", common.af=0.01){
  # Prep summary statistic collection dataframe
  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.prefix <- gsub("[ ]+", "_", paste(tolower(ref.title), metric, sep="."))

  # Overall benchmarking: all variant classes
  sub.res <- do.call("rbind", lapply(bench.dat, compress.benchmarking,
                                     min.af=common.af))
  sub.df <- data.frame(cbind(paste(ss.prefix, "common",
                                   names(bench.dat), sep="."),
                             rep(metric, 2), sub.res), row.names=NULL)
  sub.df <- rbind(sub.df,
                  c(paste(ss.prefix, "common", sep="."), metric,
                    weighted.mean(as.numeric(sub.df[, 3]),
                                  as.numeric(sub.df[, 4])),
                    sum(as.numeric(sub.df[, 4]))))
  colnames(sub.df) <- colnames(ss.df)
  ss.df <- rbind(ss.df, sub.df)

  # Benchmarking for each variant class
  vcs <- intersect(names(var.class.abbrevs),
                   unlist(lapply(bench.dat, function(l){l$df$class})))
  vc.sub.df <- do.call("rbind", lapply(vcs, function(vc){
    sub.res <- do.call("rbind", lapply(bench.dat, compress.benchmarking,
                                       vc.elig=vc, min.af=common.af))
    sub.df <- data.frame(cbind(paste(ss.prefix, "common", vc,
                                     names(bench.dat), sep="."),
                               rep(metric, 2), sub.res), row.names=NULL)
    sub.df <- rbind(sub.df,
                    c(paste(ss.prefix, "common", vc, sep="."), metric,
                      weighted.mean(as.numeric(sub.df[, 3]),
                                    as.numeric(sub.df[, 4])),
                      sum(as.numeric(sub.df[, 4]))))
  }))
  colnames(vc.sub.df) <- colnames(ss.df)
  ss.df <- rbind(ss.df, vc.sub.df)

  # Benchmarking for each variant subclass
  vscs <- intersect(names(var.subclass.abbrevs),
                    unlist(lapply(bench.dat, function(l){l$df$subclass})))
  vsc.sub.df <- do.call("rbind", lapply(vscs, function(vsc){
    sub.res <- do.call("rbind", lapply(bench.dat, compress.benchmarking,
                                       vsc.elig=vsc, min.af=common.af))
    sub.df <- data.frame(cbind(paste(ss.prefix, "common", vsc,
                                     names(bench.dat), sep="."),
                               rep(metric, 2), sub.res), row.names=NULL)
    sub.df <- rbind(sub.df,
                    c(paste(ss.prefix, "common", vsc, sep="."), metric,
                      weighted.mean(as.numeric(sub.df[, 3]),
                                    as.numeric(sub.df[, 4])),
                      sum(as.numeric(sub.df[, 4]))))
  }))
  colnames(vsc.sub.df) <- colnames(ss.df)
  ss.df <- rbind(ss.df, vsc.sub.df)

  return(ss.df)
}


######################
# Plotting Functions #
######################
# Generate a single benchmarking barplot
bench.barplot <- function(plot.dat, set.colors, title=NULL, strata.labels=NULL,
                          strata.label.srt=NULL, reverse.strata=FALSE,
                          strata.title=NULL, y.title="Variants", bar.buffer=0.2,
                          parmar=c(2, 2.5, 2, 2.5)){
  # Define plot parameters
  n.sets <- length(plot.dat)
  strata <- unique(unlist(lapply(plot.dat, function(l){rownames(l$df)})))
  n.strata <- length(strata)
  bins <- ncol(plot.dat[[1]]$df)
  xlims <- c(0, n.strata)
  ylims <- c(0, 1)
  set.counts <- do.call("rbind", lapply(strata, function(sx){
    sapply(plot.dat, function(l){if(sx %in% names(l$n)){l$n[sx]}else{0}})
  }))
  set.props <- as.data.frame(t(apply(set.counts, 1, function(v){v / sum(v)})))
  colnames(set.props) <- names(plot.dat)
  rownames(set.props) <- strata
  strata.counts <- apply(set.counts, 1, sum)
  strata.props <- strata.counts / sum(strata.counts)
  names(strata.props) <- strata
  bar.x0 <- (1:n.strata) - 1 + (bar.buffer / 2)
  bar.x1 <- (1:n.strata) - (bar.buffer / 2)
  set.pals <- lapply(set.colors, function(col.base){
    colorRampPalette(c("black", col.base, "white"))(bins + 3)[-c(1:2, bins+3)]
  })
  names(set.pals) <- names(plot.dat)

  # Set up plot area
  prep.plot.area(xlims, ylims, parmar)
  segments(x0=par("usr")[1], x1=max(bar.x1), y0=1, y1=1,
           col=annotation.color, xpd=T, lty=3)

  # Add main X axis
  if(is.null(strata.labels)){
    strata.labels <- strata
    names(strata.labels) <- strata
  }
  strata.order <- if(reverse.strata){rev(strata)}else{strata}
  sapply(strata.order, function(sx){
    si <- which(strata.order == sx)
    if(is.null(strata.label.srt)){
      axis(1, at=si-0.5, tick=F, cex.axis=4.5/6, line=-1,
           labels=parse(text=strata.labels[sx]))
    }else{
      text(x=si-0.5, y=-0.08, cex=4.5/6, srt=strata.label.srt,
           labels=parse(text=strata.labels[sx]), xpd=T)
    }
  })
  mtext(side=1, text=strata.title,
        line=if(is.null(strata.label.srt)){0.9}else{1.1})

  # Add main Y axis
  clean.axis(2, at=seq(0, 1, by=0.25), label.units="percent", label.line=-0.8)
  axis(2, at=0.5, tick=F, line=0.6, labels=y.title)
  mtext(text=title, side=3, line=0.15)

  # Add legend to right Y axis
  af.breaks <- paste("<", plot.dat[[1]]$breaks, sep="")
  right.y.cs <- cumsum(c(0, as.numeric(plot.dat[[n.sets]]$df[strata.order[n.strata], ])))
  right.y.mids <- sapply(1:length(af.breaks), function(k){mean(right.y.cs[k+(0:1)])})
  right.y.lab.pos <- yaxis.legend(af.breaks, x=max(bar.x1),
                                  y.positions=right.y.mids,
                                  sep.wex=0.05*diff(par("usr")[1:2]),
                                  upper.limit=0.5,
                                  colors=hex2grey(set.pals[[n.sets]]),
                                  label.cex=4.5/6, lwd=2, return.label.pos=T)
  axis(4, at=max(right.y.lab.pos)+0.1, las=2, tick=F, line=-0.45,
       cex.axis=5/6, labels=bquote(Delta * "AF:"))
  set.labs.y <- seq(0.95, 0.7, length.out=n.sets+1)[1:n.sets]
  points(x=rep(par("usr")[2]+(1.5*bar.buffer), 2), y=set.labs.y,
         pch=15, col=set.colors, xpd=T)
  text(x=rep(par("usr")[2]+bar.buffer, 2), y=set.labs.y-0.01,
       cex=5/6, xpd=T, labels=names(plot.dat), pos=4)

  # Add bars
  sapply(strata.order, function(sx){
    si <- which(strata.order == sx)
    x0 <- bar.x0[si]
    x1 <- bar.x1[si]
    xs <- quantile(c(x0, x1), probs=cumsum(c(0, as.numeric(set.props[sx, ]))))
    sapply(1:n.sets, function(k){
      set.name <- names(plot.dat)[k]
      set.df <- plot.dat[[k]]$df
      if(sx %in% rownames(set.df)){
        ys <- cumsum(c(0, as.numeric(set.df[sx, ])))
        rect(xleft=xs[k], xright=xs[k+1], ybottom=ys[-length(ys)], ytop=ys[-1],
             col=set.pals[[k]], border=NA, bty="n")
      }
    })
    segments(x0=x0, x1=x1, y0=0, y1=0, xpd=T, lend="butt")
  })
}

# Main wrapper function for plotting benchmarking summary results
plot.all.benchmarking <- function(bench.dat, out.prefix, set.colors,
                                  strata.labels=NULL, reverse.strata=FALSE,
                                  strata.title=NULL, ref.title="ref. data",
                                  y.title="Variants", title.prefix="Concordance of",
                                  title.preposition="with", common.af=0.01,
                                  ss.metric.name="concordance"){
  # Generate main plot with all variation
  main.dat <- lapply(bench.dat, extract.bench.res)
  main.total.n <- sum(unlist(sapply(main.dat, function(l){l$n})))
  bench.title <- paste(title.prefix, " ",
                       clean.numeric.labels(main.total.n, acceptable.decimals=1),
                       " variants\n", title.preposition, " ",
                       ref.title, sep="")
  pdf(paste(out.prefix, "all.barplot.pdf", sep="."), height=2.25, width=2.6)
  bench.barplot(main.dat, set.colors, title=bench.title,
                strata.labels=strata.labels, strata.title=strata.title,
                reverse.strata=reverse.strata, y.title=y.title)
  dev.off()

  # Generate variant class-specific subplots
  vcs <- intersect(names(var.class.abbrevs),
                   unlist(lapply(bench.dat, function(l){l$df$class})))
  for(vc in vcs){
    sub.dat <- lapply(bench.dat, extract.bench.res, vc.elig=vc)
    sub.total.n <- sum(unlist(sapply(sub.dat, function(l){l$n})))
    vc.title <- if(vc == "indel"){"indel"}else{var.class.abbrevs[vc]}
    bench.title <- paste(title.prefix, " ",
                         clean.numeric.labels(sub.total.n, acceptable.decimals=1),
                         " ", vc.title, "s\n",
                         title.preposition, " ", ref.title, sep="")
    pdf(paste(out.prefix, vc, "barplot.pdf", sep="."),
        height=2.25, width=2.6)
    bench.barplot(sub.dat, set.colors, title=bench.title,
                  strata.labels=strata.labels, strata.title=strata.title,
                  reverse.strata=reverse.strata, y.title=y.title)
    dev.off()
  }

  # Generate subclass plot for common short variants, if present
  has.short.vars <- any(unlist(lapply(bench.dat, function(l){
    sub.df <- l$df
    if(!("af" %in% colnames(sub.df))){
      return(FALSE)
    }
    any(sub.df$class %in% c("snv", "indel")
        & clean.breaks(sub.df$af) < common.af)
  })))
  if(has.short.vars){
    short.vsc.dat <- lapply(bench.dat, extract.bench.res,
                            vc.elig=c("snv", "indel"), common.af=common.af,
                            summarize.by.vsc=TRUE)
    bench.title <- paste(title.prefix, " common short\nvariants ",
                         title.preposition, " ", ref.title, sep="")
    vsc.strata.labels <- var.subclass.abbrevs[rownames(short.vsc.dat[[1]]$df)]
    pdf(paste(out.prefix, "short_subclasses.barplot.pdf", sep="."),
        height=2.25, width=2.6)
    bench.barplot(short.vsc.dat, set.colors, title=bench.title,
                  strata.labels=vsc.strata.labels,
                  strata.title="Short variant class",
                  y.title=y.title)
    dev.off()
  }

  # Generate subclass plot for common SVs, if present
  has.svs <- any(unlist(lapply(bench.dat, function(l){
    sub.df <- l$df
    if(!("af" %in% colnames(sub.df))){
      return(FALSE)
    }
    any(sub.df$class == "sv" & clean.breaks(sub.df$af) < common.af)
  })))
  if(has.svs){
    sv.vsc.dat <- lapply(bench.dat, extract.bench.res, vc.elig=c("sv"),
                         common.af=common.af, summarize.by.vsc=TRUE)
    bench.title <- paste(title.prefix, " common SVs\n",
                         title.preposition, " ", ref.title, sep="")
    all.vscs <- unique(unlist(sapply(sv.vsc.dat, function(l){rownames(l$df)})))
    vsc.strata.labels <- var.subclass.abbrevs[all.vscs]
    pdf(paste(out.prefix, "sv_subclasses.barplot.pdf", sep="."),
        height=2.25, width=2.6)
    bench.barplot(sv.vsc.dat, set.colors, title=bench.title,
                  strata.labels=vsc.strata.labels, strata.title="SV class",
                  y.title=y.title, strata.label.srt=30)
    dev.off()
  }

  # Collect & return summary metrics
  gather.bench.sumstats(bench.dat, ref.title, metric=ss.metric.name,
                        common.af=common.af)
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot site benchmarking metrics")
parser$add_argument("--sens-by-af", metavar=".tsv", action="append",
                    help=paste("Sensitivity benchmarking results by AF. Can be",
                               "specified multiple times for multiple",
                               "evaluation interval sets."))
parser$add_argument("--ppv-by-af", metavar=".tsv", action="append",
                    help=paste("PPV benchmarking results by AF. Can be",
                               "specified multiple times for multiple",
                               "evaluation interval sets."))
parser$add_argument("--set-name", metavar="string", action="append",
                    help=paste("Optional names for each set of benchmarking",
                               "results supplied"))
parser$add_argument("--ref-title", metavar="path", type="character",
                    help="String title for comparison dataset/cohort")
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
parser$add_argument("--custom-constants", metavar=".R", type="character",
                    help="Optional file of custom constants to use for plotting")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("sens_by_af" = c("~/scratch/site_bench_metric_inputs/dfci-g2c.v1.initial_qc.gnomad_v4.1.Easy.sensitivity_by_freq.merged.tsv.gz",
#                               "~/scratch/site_bench_metric_inputs/dfci-g2c.v1.initial_qc.gnomad_v4.1.Hard.sensitivity_by_freq.merged.tsv.gz"),
#              "ppv_by_af" = c("~/scratch/site_bench_metric_inputs/dfci-g2c.v1.initial_qc.gnomad_v4.1.Easy.ppv_by_freq.merged.tsv.gz",
#                               "~/scratch/site_bench_metric_inputs/dfci-g2c.v1.initial_qc.gnomad_v4.1.Hard.ppv_by_freq.merged.tsv.gz"),
#              "set_name" = c("Easy", "Hard"),
#              "ref_title" = "gnomAD v4.1",
#              "common_af" = 0.001,
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/g2c.qc.test")

# Load custom constants if optioned
if(!is.null(args$custom_constants)){
  source(args$custom_constants)
}

# Load sensitivity & PPV data
sens.by.af <- load.benchmark.tsvs(args$sens_by_af, args$set_name)
ppv.by.af <- load.benchmark.tsvs(args$ppv_by_af, args$set_name)

# Define palette for sets
set.colors <- rev(RLCtools::categorical.rainbow(2))

# Plot sensitivity vs AF
if(!is.null(sens.by.af)){
  af.strata.labels <- sapply(clean.breaks(unique(sens.by.af[[1]]$df[, 3])), function(l){
    if(log10(l) >= -2){as.character(l)}else{paste("10^", log10(l))}
  })
  strata.title <- if(is.null(args$ref_title)){
    "Allele frequency"
  }else{
    paste("AF (", args$ref_title, ")", sep="")
  }
  sens.ss <- plot.all.benchmarking(sens.by.af,
                                   paste(args$out_prefix, "sens_by_af", sep="."),
                                   set.colors, strata.labels=af.strata.labels,
                                   reverse.strata=TRUE, strata.title=strata.title,
                                   args$ref_title, y.title="Rediscovery rate",
                                   "Rediscovery of", "from", common.af=args$common_af,
                                   ss.metric.name="sensitivity")
}else{
  sens.ss <- NULL
}

# Plot PPV vs AF
if(!is.null(ppv.by.af)){
  af.strata.labels <- sapply(clean.breaks(unique(ppv.by.af[[1]]$df[, 3])), function(l){
    if(log10(l) >= -2){as.character(l)}else{paste("10^", log10(l))}
  })
  strata.title <- if(is.null(args$ref_title)){
    "Allele frequency"
  }else{
    paste("AF (this cohort)", sep="")
  }
  ppv.ss <- plot.all.benchmarking(ppv.by.af,
                                  paste(args$out_prefix, "ppv_by_af", sep="."),
                                  set.colors, strata.labels=af.strata.labels,
                                  reverse.strata=TRUE, strata.title=strata.title,
                                  args$ref_title, y.title="Confirmation rate",
                                  "Confirmation of", "versus",
                                  common.af=args$common_af, ss.metric.name="ppv")
}else{
  ppv.ss <- NULL
}

# Combine summary statistics and write to .tsv
ss.out <- do.call("rbind", list(sens.ss, ppv.ss))
for(sn in args$set_name){
  ss.out[, 1] <- gsub(paste(".", sn, "$", sep=""),
                      paste(".", tolower(sn), sep=""),
                      ss.out[, 1])
}
colnames(ss.out)[1] <- paste("#", colnames(ss.out)[1], sep="")
write.table(ss.out, paste(args$out_prefix, "summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

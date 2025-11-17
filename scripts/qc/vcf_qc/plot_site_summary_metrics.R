#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot variant site-level summary distributions for quality control


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(DescTools, quietly=TRUE)
require(G2CR, quietly=TRUE)
require(viridis, quietly=TRUE)
load.constants("all")


##################
# Data Functions #
##################
# Load all SV sizes from an SV sites BED and split by subclass
read.sv.sizes <- function(tsv.in){
  if(is.null(tsv.in)){return(NULL)}
  df <- read.table(tsv.in, header=T, sep="\t", quote="",
                   check.names=F, comment.char="")[, c("subclass", "size")]
  sv.classes.in.sites <- rev(setdiff(intersect(names(sv.colors),
                                               unique(df$subclass)),
                                     c("CTX", "BND")))
  sv.sizes <- lapply(sv.classes.in.sites, function(vsc){
    log10(as.numeric(df$size[which(df$subclass == vsc)]))
  })
  names(sv.sizes) <- sv.classes.in.sites
  return(sv.sizes)
}

# Compute allele frequency summary metrics for all classes & subclasses
get.af.ss <- function(af.d, common.af=0.01){
  af.df <- af.d$df
  breaks <- af.d$breaks

  # Get column indexes for singleton, rare non-singleton, and common
  # This assumes the first frequency column is always singleton (or close enough)
  singleton.cidx <- 3
  rare.cidx <- setdiff(which(breaks < common.af) + 2, singleton.cidx)
  common.cidx <- which(breaks >= common.af) + 2
  cidx.list <- list(singleton.cidx, rare.cidx, common.cidx)
  bin.prefixes <- c("pct_singletons", "pct_rare", "pct_common")

  # Get statistics across all variant classes
  n.all <- sapply(cidx.list, function(cidx){sum(af.df[, cidx])})
  ss.df <- data.frame("analysis" = paste(bin.prefixes, "all", sep="."),
                      "measure" = "percent",
                      "value" = n.all / sum(n.all),
                      "n" = sum(n.all))

  # Get statistics across major variant classes
  vc.ss <- do.call("rbind", lapply(names(var.class.abbrevs), function(vc){
    vc.n <- sapply(cidx.list, function(cidx){
      sum(af.df[which(af.df$class %in% vc), cidx])
    })
    if(sum(vc.n) > 0){
      data.frame("analysis" = paste(bin.prefixes, vc, sep="."),
                 "measure" = "percent",
                 "value" = vc.n / sum(vc.n),
                 "n" = sum(vc.n))
    }
  }))

  # Get statistics across variant subclasses
  vsc.ss <- do.call("rbind", lapply(names(var.subclass.abbrevs), function(vsc){
    vsc.n <- sapply(cidx.list, function(cidx){
      sum(af.df[which(af.df$subclass %in% vsc), cidx])
    })
    if(sum(vsc.n) > 0){
      data.frame("analysis" = paste(bin.prefixes, vsc, sep="."),
                 "measure" = "percent",
                 "value" = vsc.n / sum(vsc.n),
                 "n" = sum(vsc.n))
    }
  }))

  as.data.frame(rbind(ss.df, vc.ss, vsc.ss))
}


######################
# Plotting Functions #
######################
# Summary plot of variant counts by class & subclass
plot.counts.by.vsc <- function(df, has.short.variants=TRUE, has.svs=TRUE,
                               ref.size.d=NULL, ref.title=NULL,
                               bar.sep=0.1, ref.pt.cex=2/3,
                               parmar=c(0.1, 7.5, 2, 2)){

  # Get summary statistics for output .tsv
  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.df[1, ] <- c("site_count.all_variants", "count", rep(sum(df[, -(1:2)]), 2))
  vc.counts <- data.frame(do.call("rbind", lapply(names(var.class.abbrevs), function(vc){
    c(paste("site_count", vc, sep="."), "count",
      rep(sum(df[which(df$class == vc), -(1:2)]), 2))
  })))
  vsc.counts <- data.frame(do.call("rbind", lapply(names(var.subclass.abbrevs), function(vsc){
    c(paste("site_count", vsc, sep="."), "count",
      rep(sum(df[which(df$subclass == vsc), -(1:2)]), 2))
  })))
  colnames(vc.counts) <- colnames(vsc.counts) <- colnames(ss.df)
  if(has.short.variants){
    n.snv <- as.numeric(vc.counts[vc.counts$analysis == "site_count.snv", "n"])
    titv.r <- as.numeric(vsc.counts[vsc.counts$analysis == "site_count.ti", "n"]) / n.snv
    n.indel <- as.numeric(vc.counts[vc.counts$analysis == "site_count.indel", "n"])
    insdel.r <- as.numeric(vsc.counts[vsc.counts$analysis == "site_count.ins", "n"]) / n.indel
    short.ratios <- data.frame("analysis"=c("site_ratios.snv", "site_ratios.indel"),
                               "measure"=c("ti_tv_ratio", "ins_del_ratio"),
                               "value"=c(titv.r, insdel.r),
                               "n"=c(n.snv, n.indel))
  }else{
    short.ratios <- data.frame("analysis"=character(), "measure"=character(),
                               "value"=numeric(), "n"=numeric())
  }
  if(has.svs){
    n.sv <- as.numeric(vc.counts[vc.counts$analysis == "site_count.sv", "n"])
    DUPDEL.r <- sum(as.numeric(vsc.counts[vsc.counts$analysis %in% c("site_count.DUP", "site_count.INS"), "n"])) / n.sv
    sv.ratios <- data.frame("analysis"="site_ratios.sv", "measure"="DUP_DEL_ratio",
                               "value"=DUPDEL.r, "n"=n.sv)
  }else{
    sv.ratios <- data.frame("analysis"=character(), "measure"=character(),
                            "value"=numeric(), "n"=numeric())
  }
  ss.df <- as.data.frame(rbind(ss.df, vc.counts, vsc.counts, short.ratios, sv.ratios))

  # Simplify count data
  k <- log10(apply(df[, -c(1:2)], 1, sum))
  k.order <- order(k)
  k <- k[k.order]
  df <- df[k.order, ]
  names(k) <- df$subclass

  # Simplify & align reference counts, if provided
  if(is.null(ref.size.d)){
    add.ref <- FALSE
  }else{
    add.ref <- TRUE
    ref.df <- ref.size.d$df
    ref.df <- ref.df[which(ref.df$subclass %in% names(k)), ]
    ref.k <- log10(apply(ref.df[, -c(1:2)], 1, sum))
    names(ref.k) <- ref.df$subclass
    ref.k <- sapply(names(k), function(vsc){
      if(vsc %in% names(ref.k)){as.numeric(ref.k[vsc])}else{0}
    })
  }

  # Get plot parameters
  if(add.ref){
    xlims <- c(0, max(c(ceiling(k), ceiling(ref.k))))
  }else{
    xlims <- c(0, max(c(ceiling(k))))
  }
  ylims <- c(if(add.ref){-1}else{0}, length(k)+bar.sep)
  bar.cols <- var.class.colors[df$class]
  bar.labs <- sapply(10^k, function(lk){
    if(lk < 10000){
      prettyNum(lk, big.mark=",")
    }else{
      clean.numeric.labels(lk, acceptable.decimals=1, min.label.length=3)
    }
  })

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)

  # Add bars
  rect(xleft=rep(0, length(k)), xright=k,
       ybottom=(1:length(k)) - 1 + bar.sep, ytop=(1:length(k)) - bar.sep,
       col=bar.cols, xpd=T)

  # Add reference markers, if optioned
  if(add.ref){
    ref.pw.col <- remap(as.character(ref.k <= k),
                        c("TRUE"="white", "FALSE"=var.ref.color))
    segments(y0=(1:length(ref.k)) - 1 + (2.5*bar.sep),
             y1=(1:length(ref.k)) - (2.5*bar.sep),
             x0=ref.k, x1=ref.k, col=ref.pw.col, lend="round")
    ref.legend.buffer <- 0.03
    if(!is.null(ref.title)){
      ref.legend <- paste("Hashes from", ref.title)
    }else{
      ref.legend <- "Hashes from ref. cohort"
    }
    text(x=par("usr")[2]+(5*ref.legend.buffer*diff(par("usr")[1:2])),
         y=par("usr")[3]+(ref.legend.buffer*diff(par("usr")[3:4])),
         labels=ref.legend, pos=2, cex=5/6, xpd=T, col=var.ref.color, xpd=T)
  }

  # Add count labels to right Y axis
  axis(4, at=(1:length(k))-0.5, tick=F, line=-0.9,
       cex.axis=4.5/6, labels=bar.labs, las=2)

  # Add top Y axis
  minor.ticks <- log10(logscale.minor)
  minor.ticks <- minor.ticks[which(minor.ticks < max(ylims))]
  axis(3, at=minor.ticks, labels=NA, tck=-0.0125, lwd=0.5)
  major.at <- log10(logscale.major.bp)
  major.labels <- sapply(logscale.major.bp, clean.numeric.labels)
  major.labels[which(!major.labels %in% c("1", "1k", "1M", "1B"))] <- NA
  clean.axis(3, at=major.at, labels=major.labels, title="Total variants",
             label.line=-0.9, title.line=-0.1, tck=-0.03,
             infinite.positive=TRUE)

  # Add left margin labels
  axis(2, at=(1:length(k))-0.5, las=2, line=-0.9, tick=F,
       labels=var.subclass.names.short[df$subclass], cex.axis=5/6)
  vc.x <- -0.9 * diff(par("usr")[1:2])
  bracket.lab.buf <- -0.075 * vc.x
  if(has.short.variants){
    snv.bracket.y <- c(min(which(df$class == "snv")) - 1 + bar.sep,
                       max(which(df$class == "snv")) - bar.sep)
    staple.bracket(x0=vc.x, x1=vc.x, y0=snv.bracket.y[1], y1=snv.bracket.y[2],
                   staple.len=0.35, accent.len=0.2,
                   staple.accent.color=var.class.colors["snv"])
    text(x=vc.x+bracket.lab.buf, y=mean(snv.bracket.y)-0.1, labels="SNVs",
         cex=5/6, pos=2, xpd=T)

    indel.bracket.y <- c(min(which(df$class == "indel")) - 1 + bar.sep,
                         max(which(df$class == "indel")) - bar.sep)
    staple.bracket(x0=vc.x, x1=vc.x, y0=indel.bracket.y[1], y1=indel.bracket.y[2],
                   staple.len=0.35, accent.len=0.2,
                   staple.accent.color=var.class.colors["indel"])
    text(x=vc.x+bracket.lab.buf, y=mean(indel.bracket.y)-0.1, labels="Indels\n(1-49 bp)",
         cex=5/6, pos=2, xpd=T)
  }
  if(has.svs){
    sv.bracket.y <- c(min(which(df$class == "sv")) - 1 + bar.sep,
                      max(which(df$class == "sv")) - bar.sep)
    staple.bracket(x0=vc.x, x1=vc.x, y0=sv.bracket.y[1], y1=sv.bracket.y[2],
                   staple.len=0.35, accent.len=0.2,
                   staple.accent.color=var.class.colors["sv"])
    text(x=vc.x+bracket.lab.buf, y=mean(sv.bracket.y)-0.1,
         labels="Structural\nvariants\n(SVs;\n>49 bp)",
         cex=5/6, pos=2, xpd=T)
  }

  return(ss.df)
}

# Volcano of signed variant sizes
plot.size.volcano <- function(size.d, ref.size.d=NULL, ref.title=NULL,
                              snv.width=0.2, snv.gap=0.15, indel.gap=0,
                              minor.tck=-0.01, major.tck=-0.0275, ref.lty=1,
                              ref.hex.gap=0.01, parmar=c(2, 2.75, 0.25, 0.1)){
  # Get partitioned densities
  df <- size.d$df
  GAIN.vsc <- c("INS", "DUP", "CNV")
  breaks <- log10(size.d$breaks)
  snv.k <- log10(sum(as.numeric(apply(df[which(df$class == "snv"),
                                         -c(1:2)], 2, sum))))
  do.snvs <- sum(10^snv.k, na.rm=T) > 0
  ins.k <- log10(as.numeric(df[which(df$subclass == "ins"), -c(1:2)]))
  del.k <- log10(as.numeric(df[which(df$subclass == "del"), -c(1:2)]))
  do.indels <- sum(10^c(ins.k, del.k), na.rm=T) > 0
  GAIN.k <- log10(as.numeric(apply(df[which(df$subclass %in% GAIN.vsc),
                                      -c(1:2)], 2, sum)))
  LOSS.k <- log10(as.numeric(df[which(df$subclass == "DEL"), -c(1:2)]))
  do.svs <- sum(10^c(GAIN.k, LOSS.k), na.rm=T) > 0

  # Get reference densities, if optioned
  if(is.null(ref.size.d)){
    add.ref <- do.ref.snvs <- do.ref.indels <- do.ref.svs <- FALSE
  }else{
    add.ref <- TRUE
    ref.breaks <- log10(ref.size.d$breaks)
    ref.df <- ref.size.d$df
    ref.snv.k <- log10(sum(as.numeric(apply(ref.df[which(ref.df$class == "snv"),
                                                   -c(1:2)], 2, sum))))
    do.ref.snvs <- sum(10^ref.snv.k, na.rm=T) > 0
    ref.ins.k <- log10(as.numeric(ref.df[which(ref.df$subclass == "ins"), -c(1:2)]))
    ref.del.k <- log10(as.numeric(ref.df[which(ref.df$subclass == "del"), -c(1:2)]))
    do.ref.indels <- sum(10^c(ref.ins.k, ref.del.k), na.rm=T) > 0
    ref.GAIN.k <- log10(as.numeric(apply(ref.df[which(ref.df$subclass %in% GAIN.vsc),
                                                -c(1:2)], 2, sum)))
    ref.LOSS.k <- log10(as.numeric(ref.df[which(ref.df$subclass == "DEL"), -c(1:2)]))
    do.ref.svs <- sum(10^c(ref.GAIN.k, ref.LOSS.k), na.rm=T) > 0
  }

  # Get plot values
  xlims <- c(-1, 1) * (max(breaks) + 0.1 + (0.5*snv.width) + snv.gap + indel.gap)
  xlims[1] <- xlims[1]-0.1
  if(add.ref){
    ref.df.for.lims <- merge(df[, c("class", "subclass")], ref.df, all=F, sort=F)
    ylims <- c(0, 1.05*log10(max(max(df[, -c(1:2)]), max(ref.df.for.lims[, -c(1:2)]))))
  }else{
    ylims <- c(0, log10(max(df[, -c(1:2)])))
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="r")

  # Add SNVs
  if(do.snvs){
    rect(xleft=-0.5*snv.width, xright=0.5*snv.width,
         ybottom=0, ytop=snv.k, col=var.class.colors["snv"],
         border=NA, bty="n")
    if(do.ref.snvs){
        if(ref.snv.k >= (1-ref.hex.gap)*snv.k){
          snv.ref.color <-var.ref.color
        }else{
          snv.ref.color <-"white"
        }
        segments(x0=-0.5*snv.width, x1=0.5*snv.width, y0=ref.snv.k, y1=ref.snv.k,
                 lty=ref.lty, col=snv.ref.color, xpd=T)
      }
  }

  # Add indels
  if(do.indels){
    # Add indel polygons
    indel.idx <- which(breaks >= 0 & breaks < log10(50))
    ins.xy <- step.function(x=breaks[indel.idx] + snv.gap + (0.5 * snv.width),
                            y=ins.k[indel.idx], offset=1)
    polygon(x=c(ins.xy$x, rev(ins.xy$x)),
            y=c(ins.xy$y, rep(0, length(ins.xy$x))),
            col=var.class.colors["indel"], border=NA, bty="n")
    del.xy <- step.function(x=breaks[indel.idx] + snv.gap + (0.5 * snv.width),
                            y=del.k[indel.idx], offset=1)
    polygon(x=-c(del.xy$x, rev(del.xy$x)),
            y=c(del.xy$y, rep(0, length(del.xy$x))),
            col=var.class.colors["indel"], border=NA, bty="n")

    # Add indel reference ticks
    if(do.ref.indels){
        ref.indel.idx <- which(ref.breaks >= 0 & ref.breaks <= log10(50))
        ref.ins.xy <- list("x"=ref.breaks[ref.indel.idx] + snv.gap + (0.5 * snv.width),
                           "y"=ref.ins.k[ref.indel.idx])
        ref.ins.col <- sapply(1:(length(ref.ins.xy$x)-1), function(r){
          ovr.idx <- which(sapply(ins.xy$x, is.inside, interval=ref.ins.xy$x[c(r, r+1)]))
          if(length(ovr.idx) == 0){return(var.ref.color)}
          if(mean(ref.ins.xy$y[r] >= (1-ref.hex.gap)*ins.xy$y[ovr.idx]) >= 0.5){
            return(var.ref.color)
          }else{
            return("white")
          }
        })
        segments(x0=ref.ins.xy$x[-length(ref.ins.xy$x)],
                 x1=ref.ins.xy$x[-1],
                 y0=ref.ins.xy$y[-length(ref.ins.xy$y)],
                 y1=ref.ins.xy$y[-length(ref.ins.xy$y)],
                 lty=ref.lty, col=ref.ins.col, lend="butt", xpd=T)

        ref.del.xy <- list("x"=ref.breaks[ref.indel.idx] + snv.gap + (0.5 * snv.width),
                           "y"=ref.del.k[ref.indel.idx])
        ref.del.col <- sapply(1:(length(ref.del.xy$x)-1), function(r){
          ovr.idx <- which(sapply(del.xy$x, is.inside, interval=ref.del.xy$x[c(r, r+1)]))
          if(length(ovr.idx) == 0){return(var.ref.color)}
          if(mean(ref.del.xy$y[r] >= (1-ref.hex.gap)*del.xy$y[ovr.idx]) >= 0.5){
            return(var.ref.color)
          }else{
            return("white")
          }
        })
        segments(x0=-ref.del.xy$x[-length(ref.del.xy$x)],
                 x1=-ref.del.xy$x[-1],
                 y0=ref.del.xy$y[-length(ref.del.xy$y)],
                 y1=ref.del.xy$y[-length(ref.del.xy$y)],
                 lty=ref.lty, col=ref.del.col, lend="butt", xpd=T)
      }
  }

  # Add indel axes
  logscale.indels <- logscale.minor[which(logscale.minor >= 1 & logscale.minor < 50)]
  axis(1, at=log10(logscale.indels) + snv.gap + (0.5 * snv.width),
       tck=minor.tck, labels=NA, lwd=0.75)
  axis(1, at=log10(c(1, 50)) + snv.gap + (0.5 * snv.width), labels=NA, tck=0)
  axis(1, at=-log10(logscale.indels) - snv.gap - (0.5 * snv.width),
       tck=minor.tck, labels=NA, lwd=0.75)
  axis(1, at=-log10(c(1, 50)) - snv.gap - (0.5 * snv.width), labels=NA, tck=0)
  clean.axis(1, at=0:1 + snv.gap + (0.5 * snv.width), tck=major.tck,
             labels=c(1, 10), cex.axis=4.5/6, label.line=-0.9)
  clean.axis(1, at=0:-1 - snv.gap - (0.5 * snv.width), tck=major.tck,
             labels=c(1, 10), cex.axis=4.5/6, label.line=-0.9)

  # Add SVs
  if(do.svs){
    # Add SV polygons
    sv.idx <- which(breaks >= log10(50))
    gain.xy <- step.function(x=breaks[sv.idx] + snv.gap + indel.gap + (0.5*snv.width),
                             y=GAIN.k[sv.idx], offset=1)
    polygon(x=c(gain.xy$x, rev(gain.xy$x)),
            y=c(gain.xy$y, rep(0, length(gain.xy$x))),
            col=var.class.colors["sv"], border=NA, bty="n", xpd=T)
    # points(x=gain.xy$x, y=gain.xy$y, type="l", xpd=T)
    loss.xy <- step.function(x=breaks[sv.idx] + snv.gap + indel.gap + (0.5*snv.width),
                             y=LOSS.k[sv.idx], offset=1)
    polygon(x=-c(loss.xy$x, rev(loss.xy$x)),
            y=c(loss.xy$y, rep(0, length(loss.xy$x))),
            col=var.class.colors["sv"], border=NA, bty="n", xpd=T)
    # points(x=-loss.xy$x, y=loss.xy$y, type="l", xpd=T)

    # Add SV reference ticks
    if(do.ref.svs){
        ref.sv.idx <- which(ref.breaks >= log10(50))
        ref.gain.xy <- list("x"=ref.breaks[sv.idx] + snv.gap + indel.gap + (0.5*snv.width),
                            "y"=ref.GAIN.k[ref.sv.idx])
        ref.gain.col <- sapply(1:(length(ref.gain.xy$x)-1), function(r){
          ovr.idx <- which(sapply(gain.xy$x, is.inside, interval=ref.gain.xy$x[c(r, r+1)]))
          if(length(ovr.idx) == 0){return(var.ref.color)}
          if(mean(ref.gain.xy$y[r] >= (1-ref.hex.gap)*gain.xy$y[ovr.idx]) >= 0.5){
            return(var.ref.color)
          }else{
            return("white")
          }
        })
        segments(x0=ref.gain.xy$x[-length(ref.gain.xy$x)],
                 x1=ref.gain.xy$x[-1],
                 y0=ref.gain.xy$y[-length(ref.gain.xy$y)],
                 y1=ref.gain.xy$y[-length(ref.gain.xy$y)],
                 lty=ref.lty, col=ref.gain.col, lend="butt", xpd=T)

        ref.loss.xy <- list("x"=ref.breaks[sv.idx] + snv.gap + indel.gap + (0.5*snv.width),
                            "y"=ref.LOSS.k[ref.sv.idx])
        ref.loss.col <- sapply(1:(length(ref.loss.xy$x)-1), function(r){
          ovr.idx <- which(sapply(loss.xy$x, is.inside, interval=ref.loss.xy$x[c(r, r+1)]))
          if(length(ovr.idx) == 0){return(var.ref.color)}
          if(mean(ref.loss.xy$y[r] >= (1-ref.hex.gap)*loss.xy$y[ovr.idx]) >= 0.5){
            return(var.ref.color)
          }else{
            return("white")
          }
        })
        segments(x0=-ref.loss.xy$x[-length(ref.loss.xy$x)],
                 x1=-ref.loss.xy$x[-1],
                 y0=ref.loss.xy$y[-length(ref.loss.xy$y)],
                 y1=ref.loss.xy$y[-length(ref.loss.xy$y)],
                 lty=ref.lty, col=ref.loss.col, lend="butt", xpd=T)
      }
  }

  # Add SV axes
  # logscale.major.bp.labels[which(logscale.major.bp.labels == "1Mb")] <- ">1Mb"
  # logscale.major.bp.labels[which(logscale.major.bp.labels == "100kb")] <- NA
  logscale.major.bp.labels <- gsub("b|p", "", logscale.major.bp.labels)
  logscale.svs <- logscale.minor[which(logscale.minor >= 50 & logscale.minor <= 1e6)]
  sv.lbp.idx <- which(logscale.major.bp >= 50 & logscale.major.bp <= 1e6)
  sv.lab.at <- logscale.major.bp[sv.lbp.idx]
  axis(1, at=log10(logscale.svs) + snv.gap + indel.gap + (0.5*snv.width),
       tck=minor.tck, labels=NA, lwd=0.75)
  axis(1, at=-log10(logscale.svs) - snv.gap - indel.gap - (0.5*snv.width),
       tck=minor.tck, labels=NA, lwd=0.75)
  clean.axis(1, at=log10(logscale.major.bp[sv.lbp.idx]) + snv.gap + indel.gap + (0.5*snv.width),
             tck=major.tck, cex.axis=4.5/6, label.line=-0.9,
             labels=logscale.major.bp.labels[sv.lbp.idx])
  clean.axis(1, at=-log10(logscale.major.bp[sv.lbp.idx]) - snv.gap - indel.gap - (0.5*snv.width),
             tck=major.tck, cex.axis=4.5/6, label.line=-0.9,
             labels=logscale.major.bp.labels[sv.lbp.idx])

  # Add directional X axis titles
  axis(1, at=mean(c(-snv.gap - (0.5*snv.width), par("usr")[1])),
       tick=F, line=0, labels="Nucleotides lost")
  axis(1, at=mean(c(snv.gap + (0.5*snv.width), par("usr")[2])),
       tick=F, line=0, labels="Nucleotides gained")

  # Add Y axis
  axis(2, at=log10(logscale.minor), labels=NA, tck=minor.tck, lwd=0.75)
  y.ax.at <- logscale.major[which(logscale.major >= 1)]
  clean.axis(2, at=log10(y.ax.at), labels=sapply(y.ax.at, clean.numeric.labels),
             title="Variant count", title.line=0.85)

  # Add legend
  if(do.snvs){
    text(x=0, y=0.95*snv.k, labels="SNVs", pos=2,
         col=var.class.colors["snv"], xpd=T)
  }
  if(do.indels){
    text(x=-log10(2)-snv.gap-snv.width, y=max(del.xy$y),
         labels="Indels", pos=2, col=var.class.colors["indel"], xpd=T)
  }
  if(do.svs){
    text(x=-log10(100)-snv.gap-indel.gap-snv.width, y=max(loss.xy$y),
         labels="SVs", pos=2, col=var.class.colors["sv"], xpd=T)
  }
  if(add.ref){
    if(!is.null(ref.title)){
      ref.legend <- paste("Hashes from\n", ref.title)
    }else{
      ref.legend <- "Hashes from\nref. cohort"
    }
    if(do.ref.svs){
      ref.legend.y <- max(ref.gain.xy$y[which(names(ref.gain.xy$x) %in%
                                            names(which(ref.breaks >= log10(500))))],
                          na.rm=T) + (0.075*diff(par("usr")[3:4]))
    }else{
      ref.legend.y <- 0.925*par("usr")[4]
    }
    text(x=par("usr")[2], y=ref.legend.y, pos=2, cex=4.5/6,
         col=var.ref.color, xpd=T, labels=ref.legend)
  }
  # n.balcpx <- sum(df[which(df$subclass %in% c("INV", "CPX", "CTX", "OTH")), -c(1:2)])
  # balcpx.lab <- clean.numeric.labels(n.balcpx)
  # text(x=par("usr")[2], y=0.8*par("usr")[4], col=annotation.color, cex=5/6, pos=2,
  #      labels=paste("Not shown:\n", balcpx.lab, "balanced &\ncomplex SVs"), xpd=T)
}

# Ridgeplot of SV sizes
sv.size.ridgeplot <- function(sv.sizes){

  # Get summary statistics for output .tsv
  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.df[1, ] <- c("sv_size.all", "mean",
                  mean(10^unlist(sv.sizes)),
                  length(unlist(sv.sizes)))
  ss.df[2, ] <- c("sv_size.all", "median",
                  median(10^unlist(sv.sizes)),
                  length(unlist(sv.sizes)))
  vsc.ss <- data.frame(do.call("rbind", lapply(names(sv.sizes), function(vsc){
    data.frame("analysis" = paste("sv_size", vsc, sep="."),
               "measure" = c("mean", "median"),
               "value" = c(mean(10^sv.sizes[[vsc]]), median(10^sv.sizes[[vsc]])),
               "n" = length(sv.sizes[[vsc]]))
  })))
  ss.df <- as.data.frame(rbind(ss.df, vsc.ss))

  # Plot sizes
  ridgeplot(sv.sizes, breaks=seq(log10(10), log10(3000000), by=0.2),
            xlims=log10(c(10, 5000000)), x.axis.side=NA,
            names=var.subclass.abbrevs[names(sv.sizes)],
            fill=sv.colors[names(sv.sizes)],
            fancy.light.fill=adjust.color.hsb(sv.colors[names(sv.sizes)], s=-0.2, b=0.2),
            fancy.quartile.color=adjust.color.hsb(sv.colors[names(sv.sizes)], s=-0.4, b=0.4),
            border=adjust.color.hsb(sv.colors[names(sv.sizes)], b=-0.2),
            border.lwd=1.5, fancy.quartile.lwd=c(NA, 1), fancy.quartile.lend="square",
            hill.overlap=0.25, parmar=c(2, 2, 0.1, 0.1))
  axis(1, at=log10(logscale.minor), tck=-0.01, labels=NA, lwd=0.75)
  clean.axis(1, at=log10(logscale.major.bp),
             labels=logscale.major.bp.labels[seq(1, length(logscale.major.bp), 2)],
             labels.at=log10(logscale.major.bp)[seq(1, length(logscale.major.bp), 2)],
             label.line=-0.9, title.line=0, title="SV size")

  return(ss.df)
}


# Plot AF distribution per variant class
plot.af.distribs <- function(af.df, breaks, ref.af.df=NULL, colors=NULL,
                             group.names=NULL, lwd=3, y.title="Variant count",
                             common.af=0.01, parmar=c(2, 2.75, 0.25, 2.5)){
  # Prepare AF data
  af.dat <- lapply(1:nrow(af.df), function(i){
    step.function(x=breaks, y=unlist(log10(af.df[i, ])),
                  offset=0, interpolate=TRUE)
  })
  if(!is.null(ref.af.df)){
    add.ref <- TRUE
    ref.af.dat <- lapply(1:nrow(ref.af.df), function(i){
      step.function(x=breaks, y=unlist(log10(ref.af.df[i, ])),
                    offset=0, interpolate=TRUE)
    })
  }else{
    add.ref <- FALSE
  }

  # Get plot parameters
  xlims <- c(min(sapply(af.dat, function(l){l$x})) - 0.1,
             max(sapply(af.dat, function(l){l$x})))
  if(add.ref){
    ymin <- max(c(0,
                  min(c((floor(10 * min(sapply(af.dat, function(l){l$y}))) / 10) - 0.1,
                        (floor(10 * min(sapply(ref.af.dat, function(l){l$y}))) / 10) - 0.1))))
  }else{
    ymin <- max(c(0, (floor(10 * min(sapply(af.dat, function(l){l$y}))) / 10) - 0.1))
  }
  ymax <- ceiling((ceiling(10 * max(sapply(af.dat, function(l){l$y}))) / 10) + 0.25)
  ylims <- c(ymin, ymax)
  if(is.null(colors)){
    colors <- rev(greyscale.palette(length(af.dat)))
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar=parmar)

  # Add distinction between rare & common, if optioned
  if(!is.null(common.af)){
    abline(v=log10(common.af), col=annotation.color)
    text(x=log10(common.af)+(0.035*diff(par("usr")[1:2])),
         y=par("usr")[4]-(0.04*diff(par("usr")[3:4])),
         cex=5/6, pos=2, labels="Rare", xpd=T)
    text(x=log10(common.af)-(0.035*diff(par("usr")[1:2])),
         y=par("usr")[4]-(0.04*diff(par("usr")[3:4])),
         cex=5/6, pos=4, labels="Common", xpd=T)
    text(x=log10(common.af)-(0.035*diff(par("usr")[1:2])),
         y=par("usr")[4]-(0.13*diff(par("usr")[3:4])),
         labels=bquote("(AF" >= .(paste(round(100 * common.af, 1), "%)", sep=""))),
         cex=5/6, pos=4, xpd=T)
  }

  # Add ref ticks first s/t main data can be overplotted
  if(add.ref){
    sapply(1:length(ref.af.dat), function(i){
      segments(x0=ref.af.dat[[i]]$x[c(TRUE, FALSE)],
               x1=ref.af.dat[[i]]$x[c(FALSE, TRUE)],
               y0=ref.af.dat[[i]]$y[c(TRUE, FALSE)],
               y1=ref.af.dat[[i]]$y[c(TRUE, FALSE)],
               col=MixColor(colors[i], annotation.color, 0.5),
               lend="butt", lwd=1)
    })
  }

  # Add step functions and group names (if optioned)
  sapply(1:length(af.dat), function(i){
    points(af.dat[[i]]$x, af.dat[[i]]$y, col=colors[i], type="l", lwd=lwd)
  })

  # Add axes
  x.ax.at <- 0:floor(par("usr")[1])
  x.ax.labs <- head(c(1, 0.1, 0.01), length(x.ax.at))
  if(length(x.ax.at) > length(x.ax.labs)){
    x.ax.labs <- c(x.ax.labs, paste("10^", x.ax.at[-(1:length(x.ax.labs))], sep=""))
  }
  axis(1, at=log10(logscale.minor), tck=-0.01, lwd=0.75, labels=NA)
  clean.axis(1, at=x.ax.at, labels=x.ax.labs, parse.labels=TRUE, label.line=-0.9,
             title="Allele frequency (AF)", title.line=0)
  axis(2, at=log10(logscale.minor), tck=-0.01, lwd=0.75, labels=NA)
  clean.axis(2, at=log10(logscale.major.bp),
             labels=sapply(logscale.major.bp, clean.numeric.labels),
             title=y.title, title.line=0.85)

  # Add right Y-margin group labels
  if(!is.null(group.names)){
    yaxis.legend(group.names, x=par("usr")[2]-(0.05*diff(par("usr")[1:2])),
                 y.positions=sapply(af.dat, function(l){tail(l$y, 1)}),
                 sep.wex=0, min.label.spacing=0.1*diff(par("usr")[3:4]),
                 lower.limit=par("usr")[3]+0.025*(diff(par("usr")[3:4])),
                 upper.limit=par("usr")[4]-0.025*(diff(par("usr")[3:4])),
                 colors=NA, label.colors=colors)
  }
}

# Stacked barplot of AF bins by variant size
plot.size.by.af <- function(joint.d, bar.sep=0.1, parmar=c(2.6, 2.75, 0.25, 3.75)){
  # Convert data into proportions
  df <- joint.d$df
  af.bins <- joint.d$breaks
  size.bins <- sort(unique(clean.breaks(df$size)))
  names(size.bins) <- paste("ge", size.bins, sep="")
  plot.df <- do.call("cbind", lapply(size.bins, function(ge){
    k.by.af <- apply(df[which(df$size == paste("ge", ge, sep="")), -c(1:3)], 2, sum)
    k.by.af / sum(k.by.af)
  }))

  # Get plot parameters
  xlims <- c(0, length(size.bins))
  ylims <- c(0, 1)
  af.pal <- rev(viridis(length(af.bins), begin=0.05, end=0.95))

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="r")

  # Add bottom bin labels
  x.lab.menu <- c("SNV", logscale.major.bp.labels)
  names(x.lab.menu) <- paste("ge", c(0, logscale.major.bp), sep="")
  x.labs <- x.lab.menu[names(size.bins)]
  sapply(1:length(size.bins), function(x){
    text(x=x+0.1, y=-0.035, srt=40, pos=2, labels=x.labs[x], xpd=T, cex=5/6)
  })
  mtext(1, line=1.75, text="Variant size")

  # Add left Y axis
  clean.axis(2, at=seq(0, 1, 0.25), label.units="percent",
             title="Proportion of variants", title.line=0.6, label.line=-0.75)

  # Add right Y axis/legend hybrid
  af.labels.pct <- c("\"AF\">\"10%\"", "\"AF\"<\"10%\"", "\"AF\"<\"1%\"")
  if(length(af.bins) > length(af.labels.pct)){
    af.labels <- rev(c(af.labels.pct,
                       paste("AF<10^", -(length(af.labels.pct):(length(af.bins)-1)), sep="")))
    parse.any <- TRUE
  }else{
    af.labels <- head(af.labels.pct, length(af.bins))
    parse.any <- FALSE
  }
  legend.y.pos <- (c(0, cumsum(plot.df[-nrow(plot.df), ncol(plot.df)]))
                   + cumsum(plot.df[, ncol(plot.df)])) / 2
  legend.y.pos <- yaxis.legend(af.labels, x=ncol(plot.df)-bar.sep,
                               y.positions=legend.y.pos, parse.labels=parse.any,
                               upper.limit=0.975, colors=af.pal, sep.wex=0.5,
                               label.cex=5/6, min.label.spacing=0.1,
                               return.label.pos=TRUE)

  # Add bars last
  sapply(1:length(size.bins), function(x){
    rect(xleft=rep(x-1+bar.sep, nrow(plot.df)),
         xright=rep(x-bar.sep, nrow(plot.df)),
         ybottom=c(0, cumsum(plot.df[-nrow(plot.df), x])),
         ytop=cumsum(plot.df[, x]),
         col=af.pal, border=af.pal, lwd=0.5)
    rect(xleft=x-1+bar.sep, xright=x-bar.sep, ybottom=0, ytop=1, col=NA, xpd=T)
  })
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot variant site metric distributions")
parser$add_argument("--size-distrib", metavar=".tsv", type="character",
                    help="Precomputed binned variant size distribution")
parser$add_argument("--af-distrib", metavar=".tsv", type="character",
                    help="Precomputed binned variant AF distribution")
parser$add_argument("--joint-distrib", metavar=".tsv", type="character",
                    help="Precomputed variant size vs. AF distribution")
parser$add_argument("--ref-size-distrib", metavar=".tsv", type="character",
                    help=paste("Precomputed binned variant size distribution ",
                               "for a desired external reference dataset"))
parser$add_argument("--ref-af-distrib", metavar=".tsv", type="character",
                    help=paste("Precomputed binned variant AF distribution ",
                               "for a desired external reference dataset"))
parser$add_argument("--ref-title", metavar="string", type="character",
                    help="String title for --ref-size-distrib / --ref-af-distrib")
parser$add_argument("--sv-sites", metavar=".bed", type="character",
                    help="SV sites .bed file")
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
parser$add_argument("--custom-constants", metavar=".R", type="character",
                    help="Optional file of custom constants to use for plotting")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("size_distrib" = "~/scratch/site_bench_inputs_v2_may14/dfci-g2c.v1.initial_qc.size_distribution.merged.tsv.gz",
#              "af_distrib" = "~/scratch/site_bench_inputs_v2_may14/dfci-g2c.v1.initial_qc.af_distribution.merged.tsv.gz",
#              "joint_distrib" = "~/scratch/site_bench_inputs_v2_may14/dfci-g2c.v1.initial_qc.size_vs_af_distribution.merged.tsv.gz",
#              "sv_sites" = "~/scratch/site_bench_inputs_v2_may14/dfci-g2c.v1.initial_qc.all_svs.bed.gz",
#              "common_af" = 0.001,
#              "ref_size_distrib" = "~/scratch/site_bench_inputs_v2_may14/gnomAD_v4.1.size_distribution.merged.tsv.gz",
#              "ref_af_distrib" = "~/scratch/site_bench_inputs_v2_may14/gnomAD_v4.1.af_distribution.merged.tsv.gz",
#              "ref_title" = "gnomAD v4.1",
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/g2c.qc.test")

# # YL DEV:
# args <- list("size_distrib" = "~/scratch/YL.v1.2.size_distribution.merged.tsv.gz",
#              "af_distrib" = "~/scratch/YL.v1.2.af_distribution.merged.tsv.gz",
#              "joint_distrib" = "~/scratch/YL.v1.2.size_vs_af_distribution.merged.tsv.gz",
#              "sv_sites" = "~/scratch/YL.v1.2.all_svs.bed.gz",
#              "common_af" = 0.01,
#              "ref_size_distrib" = "~/scratch/site_bench_inputs_v2_may14/gnomAD_v4.1.size_distribution.merged.tsv.gz",
#              "ref_af_distrib" = "~/scratch/site_bench_inputs_v2_may14/gnomAD_v4.1.af_distribution.merged.tsv.gz",
#              "ref_title" = "gnomAD v4.1",
#              "custom_constants" = "~/Desktop/Collins/VanAllen/jackie_younglung/data/yl_qc.custom_constants.R",
#              "out_prefix" = "~/scratch/yl.qc.test")

# Load custom constants if optioned
if(!is.null(args$custom_constants)){
  source(args$custom_constants)
}

# Read distributions
size.d <- read.compressed.distrib(args$size_distrib)
af.d <- read.compressed.distrib(args$af_distrib)
joint.d <- read.compressed.distrib(args$joint_distrib, key.cols=1:3)
sv.sizes <- read.sv.sizes(args$sv_sites)
ref.size.d <- read.compressed.distrib(args$ref_size_distrib)
ref.af.d <- read.compressed.distrib(args$ref_af_distrib)

# Check if short variant & SV data are present in compressed distributions
has.short.variants <- (any(c("snv", "indel") %in% size.d$df$class)
                       | any(c("snv", "indel") %in% af.d$df$class)
                       | any(c("snv", "indel") %in% joint.d$df$class))
has.svs <- ("sv" %in% size.d$df$class
            | "sv" %in% af.d$df$class
            | "sv" %in% joint.d$df$class)


# Summary of variant counts
count.use.d <- head(c(size.d, af.d), 1)
count.use.ref.d <- head(c(ref.size.d, ref.af.d), 1)
if(!is.null(count.use.d)){
  for(has.ref in unique(c(FALSE, !is.null(count.use.ref.d)))){
    if(has.ref){
      pdf.suffix <- "with_ref.pdf"
      ref.d.sub <- count.use.ref.d
    }else{
      pdf.suffix <- "pdf"
      ref.d.sub <- NULL
    }
    pdf(paste(args$out_prefix, "variant_count_bars", pdf.suffix, sep="."),
        height=2.25, width=2.85)
    count.ss <- plot.counts.by.vsc(count.use.d$df, has.short.variants, has.svs,
                                   ref.d.sub, ref.title=args$ref_title)
    dev.off()
  }
}else{
  count.ss <- NULL
}


# Volcano of variant sizes
if(!is.null(size.d)){
  for(has.ref in unique(c(FALSE, !is.null(ref.size.d)))){
    if(has.ref){
      pdf.suffix <- "with_ref.pdf"
      ref.d.sub <- ref.size.d
    }else{
      pdf.suffix <- "pdf"
      ref.d.sub <- NULL
    }
    pdf(paste(args$out_prefix, "variant_size_volcano", pdf.suffix, sep="."),
        height=2.25, width=4.1)
    plot.size.volcano(size.d, ref.d.sub, ref.title=args$ref_title)
    dev.off()
  }
}
# Ridgeplot of SV sizes by subclass
if(!is.null(sv.sizes)){
  pdf(paste(args$out_prefix, "sv_size_ridgeplot.pdf", sep="."),
      height=2.25, width=2.25)
  sv.size.ss <- sv.size.ridgeplot(sv.sizes)
  dev.off()
}else{
  sv.size.ss <- NULL
}


# AF plots
if(!is.null(af.d)){
  # Gather AF summary statistics
  af.ss <- get.af.ss(af.d)

  # Generate AF distribution plots
  for(has.ref in unique(c(FALSE, !is.null(ref.af.d)))){
    if(has.ref){
      pdf.suffix <- "with_ref.pdf"
    }else{
      pdf.suffix <- "pdf"
    }

    # Step function of AFs
    classes.in.af <- intersect(c("snv", "indel", "sv"), af.d$df$class)
    main.af.df <- do.call("rbind", lapply(classes.in.af, function(vc){
      apply(af.d$df[which(af.d$df$class == vc), -c(1:2)], 2, sum)
    }))
    if(has.ref){
      main.ref.af.df <- do.call("rbind", lapply(classes.in.af, function(vc){
        apply(ref.af.d$df[which(ref.af.d$df$class == vc), -c(1:2)], 2, sum)
      }))
    }else{
      main.ref.af.df <- NULL
    }
    pdf(paste(args$out_prefix, "af_distribs", pdf.suffix, sep="."),
        height=2.25, width=2.5)
    plot.af.distribs(main.af.df, breaks=log10(af.d$breaks), main.ref.af.df,
                     colors=var.class.colors[classes.in.af],
                     group.names=paste(var.class.abbrevs[classes.in.af], "s", sep=""),
                     common.af=args$common_af)
    dev.off()

    # Supplementary step function of AFs for short variants by subclass
    if(has.short.variants){
      short.af.df <- do.call("rbind",
                             lapply(c("ti", "tv", "del", "ins"), function(vsc){
                               apply(af.d$df[which(af.d$df$subclass == vsc),
                                             -c(1:2)], 2, sum)
                             }))
      if(has.ref){
        short.ref.af.df <- do.call("rbind",
                                   lapply(c("ti", "tv", "del", "ins"), function(vsc){
                                     apply(ref.af.d$df[which(ref.af.d$df$subclass == vsc),
                                                       -c(1:2)], 2, sum)
                                   }))
      }else{
        short.ref.af.df <- NULL
      }
      pdf(paste(args$out_prefix, "af_distribs.short_vars", pdf.suffix, sep="."),
          height=2.25, width=2.5)
      plot.af.distribs(short.af.df, breaks=log10(af.d$breaks), short.ref.af.df,
                       colors=var.subclass.colors[c("ti", "tv", "del", "ins")],
                       group.names=var.subclass.abbrevs[c("ti", "tv", "del", "ins")],
                       y.title="Short variant count", lwd=2, common.af=args$common_af)
      dev.off()
    }

    # Supplementary step function of AFs for SVs by subclass
    if(has.svs){
      svs.in.af <- intersect(names(sv.colors), af.d$df$subclass)
      sv.af.df <- do.call("rbind", lapply(svs.in.af, function(vsc){
        apply(af.d$df[which(af.d$df$subclass == vsc), -c(1:2)], 2, sum)
      }))
      if(has.ref){
        sv.ref.af.df <- do.call("rbind", lapply(svs.in.af, function(vsc){
          apply(ref.af.d$df[which(ref.af.d$df$subclass == vsc), -c(1:2)], 2, sum)
        }))
      }else{
        sv.ref.af.df <- NULL
      }
      pdf(paste(args$out_prefix, "af_distribs.svs", pdf.suffix, sep="."),
          height=2.25, width=2.5)
      plot.af.distribs(sv.af.df, breaks=log10(af.d$breaks), sv.ref.af.df,
                       colors=sv.colors[svs.in.af],
                       group.names=var.subclass.abbrevs[svs.in.af],
                       y.title="SV count", lwd=2, common.af=args$common_af,
                       parmar=c(2, 2.75, 0.25, 2.75))
      dev.off()
    }
  }
}else{
  af.ss <- NULL
}


# 2D comparison of sizes vs. AFs
if(!is.null(joint.d)){
  pdf(paste(args$out_prefix, "size_vs_af_bars.pdf", sep="."),
      height=2.25, width=3.25)
  plot.size.by.af(joint.d)
  dev.off()
}


# Combine summary statistics and write to .tsv
ss.out <- do.call("rbind", list(count.ss, sv.size.ss, af.ss))
colnames(ss.out)[1] <- paste("#", colnames(ss.out)[1], sep="")
write.table(ss.out, paste(args$out_prefix, "summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)


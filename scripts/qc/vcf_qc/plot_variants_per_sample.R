#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot distribution of variants per genome


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
require(zoo, quietly=TRUE)
load.constants("all")


##################
# Data Functions #
##################
# Load and process genotype counts per sample
load.gt.counts <- function(gt.dist.in){
  # Read data
  gt <- read.table(gt.dist.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(gt) <- gsub("#", "", colnames(gt))

  # For purposes of vizualization, lump "other" GTs in with homozygotes
  gt$hom <- gt$hom + gt$other
  gt$other <- NULL

  # Add collapsed parent categories per variant class
  vcs <- unique(gt$class)
  sids  <- unique(gt$sample)
  add.rows <- as.data.frame(do.call("rbind", lapply(vcs, function(vc){
    do.call("rbind", lapply(sids, function(sid){
      sub.df <- gt[which(gt$sample == sid & gt$class == vc), ]
      c(sid, vc, NA, NA, apply(sub.df[, c("het", "hom")], 2, sum))
    }))
  })))
  colnames(add.rows) <- colnames(gt)
  gt <- as.data.frame(rbind(gt, add.rows))
  gt[, c("het", "hom")] <- apply(gt[, c("het", "hom")], 2, as.numeric)

  # Add parent row for total variation per genome across all classes
  total.rows <- as.data.frame(do.call("rbind", lapply(sids, function(sid){
    sub.df <- gt[which(gt$sample == sid & gt$class %in% vcs & is.na(gt$subclass)), ]
    c(sid, "all", NA, NA, apply(sub.df[, c("het", "hom")], 2, sum))
  })))
  colnames(total.rows) <- colnames(gt)
  gt <- as.data.frame(rbind(gt, total.rows))
  gt[, c("het", "hom")] <- apply(gt[, c("het", "hom")], 2, as.numeric)

  return(gt)
}

# Load sample labels and subset to a list of eligible samples
load.labels <- function(tsv.in, elig.samples=NULL){
  if(is.null(tsv.in)){
    return(NULL)
  }
  labs <- read.table(tsv.in, header=F, sep="\t")
  if(!is.null(elig.samples)) {
    labs <- labs[which(labs[, 1] %in% elig.samples), ]
  }
  l.n <- as.character(labs[, 1])
  l.v <- as.character(labs[, 2])
  names(l.v) <- l.n

  return(l.v)
}

# Filter GT counts based on column values
filter.gt.counts <- function(gt.counts, vc=NULL, vsc=NULL, freq=NULL){
  # Subset by class
  if(!is.null(vc)){
    gt.counts <- gt.counts[which(gt.counts$class %in% vc), ]
  }

  # Subset by subclass
  if(!is.null(vsc)){
    if(is.na(vsc)){
      gt.counts <- gt.counts[which(is.na(gt.counts$subclass)), ]
    }else{
      gt.counts <- gt.counts[which(gt.counts$subclass %in% vsc), ]
    }
  }

  # Subset by frequency
  if(!is.null(freq)){
    gt.counts <- gt.counts[which(gt.counts$freq_bin %in% freq), ]
  }

  gt.counts
}

# Compute a summary metric across a subset of gt counts
calc.gt.ss <- function(gt.counts, samples, zyg=NA, summary.fx=median,
                       calc.heterozygosity=FALSE){
  # Calculate heterozygosity, if needed
  if(calc.heterozygosity){
    gt.counts$hzg <- gt.counts$het / apply(gt.counts[, c("het", "hom")], 1, sum)
  }
  # Sum counts (or calc heterozygosity) for each sample
  k <- sapply(samples, function(sid){
    s.k <- gt.counts[which(gt.counts$sample == sid),
                     if(calc.heterozygosity){"hzg"}else
                       if(is.na(zyg)){c("het", "hom")}else{zyg}]
    if(length(s.k) == 0){0}else{sum(as.numeric(unlist(s.k)))}
  })
  if(is.null(summary.fx)){
    return(k)
  }else{
    summary.fx(k)
  }
}

# Main wrapper function to collect all per-sample summary stats
gather.count.sumstats <- function(gt.counts, samples, pop=NULL, pheno=NULL){
  ss.df <- as.data.frame(do.call("rbind", lapply(unique(gt.counts$class), function(vc){
    vc.prefix <- paste("variants_per_genome", vc, sep=".")
    vc.df <- filter.gt.counts(gt.counts, vc=vc)

    do.call("rbind", lapply(unique(vc.df$subclass), function(vsc){
      vsc.prefix <- gsub("\\.NA$", "", paste(vc.prefix, vsc, sep="."))
      vsc.df <- filter.gt.counts(vc.df, vsc=vsc)

      do.call("rbind", lapply(unique(vsc.df$freq_bin), function(freq){
        freq.prefix <- gsub("\\.AF_NA$", "",
                            paste(vsc.prefix, paste("AF", freq, sep="_"), sep="."))
        freq.df <- filter.gt.counts(vsc.df, freq=freq)

        do.call("rbind", lapply(c(NA, "het", "hom"), function(zyg){
          zyg.prefix <- sub("\\.NA$", "", paste(freq.prefix, zyg, sep="."))
          med.v <- calc.gt.ss(freq.df, samples=samples, zyg=zyg, summary.fx=median)
          phi.v <- calc.gt.ss(freq.df, samples=samples, zyg=zyg,
                              summary.fx=RLCtools::dynamic.range)
          zyg.sub.res <- data.frame("analysis"=zyg.prefix,
                                    "measure"=c("median", "dynamic_range"),
                                    "value"=c(med.v, phi.v), "n"=length(samples))
          if(is.na(zyg)){
            med.hzg <- calc.gt.ss(freq.df, samples=samples, summary.fx=median,
                                calc.heterozygosity=TRUE)
            phi.hzg <- calc.gt.ss(freq.df, samples=samples,
                                  summary.fx=RLCtools::dynamic.range,
                                  calc.heterozygosity=TRUE)
            zyg.hzg.res <- data.frame("analysis"=gsub("variants_per_genome", "heterozygosity", zyg.prefix),
                                      "measure"=c("median", "dynamic_range"),
                                      "value"=c(med.hzg, phi.hzg),
                                      "n"=length(samples))
            zyg.sub.res <- rbind(zyg.sub.res, zyg.hzg.res)
          }

          if(!is.null(pop)){
            zyg.sub.add <- do.call("rbind", lapply(unique(pop), function(p){
              pop.prefix <- gsub("\\.NA$", "", paste(zyg.prefix, p, sep="."))
              pop.samples <- names(pop)[which(pop == p)]
              med.v <- calc.gt.ss(freq.df, samples=pop.samples, zyg=zyg, summary.fx=median)
              phi.v <- calc.gt.ss(freq.df, samples=pop.samples, zyg=zyg,
                                  summary.fx=RLCtools::dynamic.range)
              zyg.sub.add <- data.frame("analysis"=pop.prefix,
                                        "measure"=c("median", "dynamic_range"),
                                        "value"=c(med.v, phi.v),
                                        "n"=length(pop.samples))
              if(is.na(zyg)){
                med.hzg <- calc.gt.ss(freq.df, samples=pop.samples, summary.fx=median,
                                      calc.heterozygosity=TRUE)
                phi.hzg <- calc.gt.ss(freq.df, samples=pop.samples,
                                      summary.fx=RLCtools::dynamic.range,
                                      calc.heterozygosity=TRUE)
                zyg.hzg.res <- data.frame("analysis"=gsub("variants_per_genome", "heterozygosity", pop.prefix),
                                          "measure"=c("median", "dynamic_range"),
                                          "value"=c(med.hzg, phi.hzg),
                                          "n"=length(pop.samples))
                zyg.sub.add <- rbind(zyg.sub.add, zyg.hzg.res)
              }

              if(!is.null(pheno)){
                zyg.sub.add.pw <- do.call("rbind", lapply(unique(pheno), function(ph){
                  pw.prefix <- gsub("\\.NA$", "", paste(pop.prefix, ph, sep="."))
                  pheno.samples <- names(pheno)[which(pheno == ph)]
                  pw.samples <- intersect(pop.samples, pheno.samples)
                  if(length(pw.samples) == 0){
                    return(NULL)
                  }
                  med.v <- calc.gt.ss(freq.df, samples=pw.samples, zyg=zyg, summary.fx=median)
                  phi.v <- calc.gt.ss(freq.df, samples=pw.samples, zyg=zyg,
                                      summary.fx=RLCtools::dynamic.range)
                  zyg.sub.add.pw <- data.frame("analysis"=pw.prefix,
                                               "measure"=c("median", "dynamic_range"),
                                               "value"=c(med.v, phi.v),
                                               "n"=length(pw.samples))
                  if(is.na(zyg)){
                    med.hzg <- calc.gt.ss(freq.df, samples=pw.samples, summary.fx=median,
                                          calc.heterozygosity=TRUE)
                    phi.hzg <- calc.gt.ss(freq.df, samples=pw.samples,
                                          summary.fx=RLCtools::dynamic.range,
                                          calc.heterozygosity=TRUE)
                    zyg.hzg.res <- data.frame("analysis"=gsub("variants_per_genome", "heterozygosity", pw.prefix),
                                              "measure"=c("median", "dynamic_range"),
                                              "value"=c(med.hzg, phi.hzg),
                                              "n"=length(pw.samples))
                    zyg.sub.add.pw <- rbind(zyg.sub.add.pw, zyg.hzg.res)
                  }
                }))
                zyg.sub.add <- rbind(zyg.sub.add, zyg.sub.add.pw)
              }
              return(zyg.sub.add)
            }))
            zyg.sub.res <- rbind(zyg.sub.res, zyg.sub.add)
          }

          if(!is.null(pheno)){
            zyg.sub.add <- do.call("rbind", lapply(unique(pheno), function(p){
              pheno.prefix <- gsub("\\.NA$", "", paste(zyg.prefix, p, sep="."))
              pheno.samples <- names(pheno)[which(pheno == p)]
              med.v <- calc.gt.ss(freq.df, samples=pheno.samples, zyg=zyg, summary.fx=median)
              phi.v <- calc.gt.ss(freq.df, samples=pheno.samples, zyg=zyg,
                                  summary.fx=RLCtools::dynamic.range)
              zyg.sub.add.ph <- data.frame("analysis"=pheno.prefix,
                                           "measure"=c("median", "dynamic_range"),
                                           "value"=c(med.v, phi.v),
                                           "n"=length(pheno.samples))
              if(is.na(zyg)){
                med.hzg <- calc.gt.ss(freq.df, samples=pheno.samples, summary.fx=median,
                                      calc.heterozygosity=TRUE)
                phi.hzg <- calc.gt.ss(freq.df, samples=pheno.samples,
                                      summary.fx=RLCtools::dynamic.range,
                                      calc.heterozygosity=TRUE)
                zyg.hzg.res <- data.frame("analysis"=gsub("variants_per_genome", "heterozygosity", pheno.prefix),
                                          "measure"=c("median", "dynamic_range"),
                                          "value"=c(med.hzg, phi.hzg),
                                          "n"=length(pheno.samples))
                zyg.sub.add.ph <- rbind(zyg.sub.add.ph, zyg.hzg.res)
              }
              return(zyg.sub.add.ph)
            }))
            zyg.sub.res <- rbind(zyg.sub.res, zyg.sub.add)
          }
          return(zyg.sub.res)
        }))

      }))
    }))
  })))
  return(unique(ss.df))
}


######################
# Plotting Functions #
######################
# Helper function to plot a single count waterfall; called within main.waterfall()
plot.count.waterfall <- function(gt.counts, vc, pop=NULL, pheno=NULL,
                                 samples.sorted=NULL, pop.space.wex=0.025,
                                 pheno.space.wex=0.005, global.scaling.cex=1,
                                 parmar=c(0.25, 3.5, 0.5, 2.5)){
  # Format plotting data
  counts <- gt.counts[which(gt.counts$class == vc
                            & is.na(gt.counts$subclass)),
                      c("sample", "hom", "het")]
  rownames(counts) <- counts$sample
  if(!is.null(samples.sorted)){
    counts <- counts[samples.sorted, ]
  }
  counts$sample <- NULL
  s.feats <- data.frame("order" = 1:nrow(counts), row.names=rownames(counts))
  s.feats$pop <- if(is.null(pop)){NA}else{pop[rownames(counts)]}
  s.feats$pheno <- if(is.null(pop)){NA}else{pheno[rownames(counts)]}
  pop.breaks.x <- which(sapply(2:nrow(s.feats), function(r){
    s.feats$pop[r] != s.feats$pop[r-1]
  }))
  pheno.breaks.x <- which(sapply(2:nrow(s.feats), function(r){
    (s.feats$pop[r] == s.feats$pop[r-1]
     & s.feats$pheno[r] != s.feats$pheno[r-1])
  }))
  n.samples <- length(rownames(counts))
  total.per.sample <- apply(counts, 1, sum, na.rm=T)

  # Get plot parameters
  pop.spacer <- n.samples * pop.space.wex
  pheno.spacer <- n.samples * pheno.space.wex
  xlims <- c(-0.01*n.samples,
             n.samples
             + (pop.spacer * length(pop.breaks.x))
             + (pheno.spacer * length(pheno.breaks.x)))
  ylims <- c(0, max(total.per.sample))
  hom.col <- adjust.color.hsb(var.class.colors[vc], s=0.01, b=-0.03)
  het.col <- adjust.color.hsb(var.class.colors[vc], s=-0.01, b=0.06)

  # Determine x position for each sample
  prev.x <- 0
  sample.xleft <- c()
  for(i in 1:n.samples){
    if(i %in% pop.breaks.x){
      prev.x <- prev.x + pop.spacer
    }else if(i %in% pheno.breaks.x){
      prev.x <- prev.x + pheno.spacer
    }
    sample.xleft[samples.sorted[i]] <- prev.x
    prev.x <- prev.x + 1
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=annotation.color)

  # Add bars
  sapply(samples.sorted, function(sid){
    rect(xleft=sample.xleft[sid],
         xright=sample.xleft[sid]+1,
         ybottom=c(0, counts[sid, 1]),
         ytop=cumsum(as.numeric(counts[sid, 1:2])),
         col=c(hom.col, het.col), border=c(hom.col, het.col), lwd=0.25)
  })

  # Add groupwise medians, divided by population if provided
  sapply(unique(s.feats$pop), function(pop){
    sids <- rownames(s.feats)[which(s.feats$pop == pop)]
    group.df <- data.frame("x" = sample.xleft[sids]+0.5,
                           "k" = total.per.sample[sids],
                           row.names=sids)
    group.df <- group.df[order(group.df$x), ]
    lab.x.at <- median(group.df$x)
    # lab.y.at <- max(c(0.55*diff(par("usr")[3:4]),
    #                   quantile(group.df$k, prob=0.99)))
    lab.y.at <- median(group.df$k, na.rm=T)
    group.k.med <- median(group.df$k, na.rm=T)
    if(group.k.med < 10000){
      med.lab <- prettyNum(round(group.k.med, 0), big.mark=",")
    }else{
      med.lab <- clean.numeric.labels(round(group.k.med, 0),
                                      acceptable.decimals=2, min.label.length=2)
    }
    lab.width <- 1.075*strwidth(med.lab, cex=global.scaling.cex*5/6)
    lab.height <- 1.2*strheight(med.lab, cex=global.scaling.cex*5/6)
    rect(xleft=lab.x.at-(lab.width/2), xright=lab.x.at+(lab.width/2),
         ybottom=lab.y.at-(lab.height/2), ytop=lab.y.at+(lab.height/2),
         col=adjustcolor("white", alpha=3/4), border=NA)
    # text(x=lab.x.at, y=lab.y.at, cex=global.scaling.cex*5/6,
    #      xpd=T, col="white", font=2, labels=med.lab)
    text(x=lab.x.at, y=lab.y.at, cex=global.scaling.cex*5/6, xpd=T, labels=med.lab,
         col=if(is.na(pop)){"black"}else{pop.palettes[[pop]]["dark1"]})
  })

  # Add left Y axis
  clean.axis(2, title=paste(var.class.abbrevs[vc], "s", sep=""),
             label.units="count", infinite.positive=T, min.ticks=3, max.ticks=4,
             title.line=parmar[2]-2.1, cex.title=global.scaling.cex*7.5/6,
             cex.axis=global.scaling.cex)

  # Add right Y axis-legend hybrid
  legend.y.at <- apply(sapply(tail(ceiling(0.01*n.samples):n.samples, 25), function(si){
    rollmean(cumsum(c(0, as.numeric(counts[samples.sorted[si], 1:2]))), k=2)
  }), 1, mean)
  yaxis.legend(c("Hom.", "Het."), x=max(sample.xleft)+1, y.positions=legend.y.at,
               sep.wex=0.005*diff(par("usr")[1:2]),
               min.label.spacing=0.075*diff(par("usr")[3:4]),
               colors=c(hom.col, het.col), label.cex=5/6*global.scaling.cex)
}

# Helper function to add groupwise markers to waterfall plot
add.waterfall.markers <- function(order.df, pop.map, pheno.map,
                                  pop.space.wex=0.025, pheno.space.wex=0.005,
                                  group.labels=c("Ancestry", "Pheno."),
                                  rect.hex.buffer=0.125, global.scaling.cex=1,
                                  parmar=c(1.35, 3.5, 0.25, 2.5)){
  # Determine breaks between groups
  s.feats <- data.frame("order" = 1:nrow(order.df), row.names=rownames(order.df))
  s.feats$pop <- if(is.null(pop)){NA}else{pop[rownames(order.df)]}
  s.feats$pheno <- if(is.null(pop)){NA}else{pheno[rownames(order.df)]}
  n.samples <- nrow(s.feats)
  pop.spacer <- pop.space.wex * n.samples
  pheno.spacer <- pheno.space.wex * n.samples
  last.x <- 0
  pop.labs <- pheno.labs <- pop.starts <- pop.ends <- pheno.starts <- pheno.ends <- c()
  for(r in 2:n.samples){
    last.x <- last.x + 1
    if(s.feats$pop[r] != s.feats$pop[r-1]){
      pop.ends <- c(pop.ends, last.x)
      pheno.ends <- c(pheno.ends, last.x)
      last.x <- last.x + pop.spacer
      pop.starts <- c(pop.starts, last.x)
      pheno.starts <- c(pheno.starts, last.x)
      pop.labs <- c(pop.labs, s.feats$pop[r-1])
      pheno.labs <- c(pheno.labs, s.feats$pheno[r-1])
    }else{
      if(s.feats$pheno[r] != s.feats$pheno[r-1]){
        pheno.ends <- c(pheno.ends, last.x)
        last.x <- last.x + pheno.spacer
        pheno.starts <- c(pheno.starts, last.x)
        pheno.labs <- c(pheno.labs, s.feats$pheno[r-1])
      }
    }
    if(r == n.samples){
      if(length(pop.ends) > 0){
        pop.ends <- c(pop.ends, last.x)
        pop.labs <- c(pop.labs, s.feats$pop[r])
      }
      if(length(pheno.ends) > 0){
        pheno.ends <- c(pheno.ends, last.x)
        pheno.labs <- c(pheno.labs, s.feats$pheno[r])
      }
    }
  }
  do.pop <- length(pop.starts) > 0
  if(do.pop){
    pop.starts <- c(0, pop.starts)
  }
  do.pheno <- length(pheno.starts) > 0
  if(do.pheno){
    pheno.starts <- c(0, pheno.starts)
  }

  # Get plot parameters
  xlims <- c(-0.01*n.samples, max(c(pop.ends, pheno.ends, n.samples)))
  ylims <- c(0, ncol(order.df)-1)
  if(all(pop.labs %in% names(pop.colors))){
    pop.pal <- pop.colors[names(pop.map)]
    pop.names <- pop.names.short[names(pop.map)]
  }else{
    pop.pal <- categorical.rainbow(length(pop.map))
    pop.names <- names(pop.map)
    names(pop.names) <- names(pop.pal) <- names(pop.map)
  }
  if(all(pheno.labs %in% names(cancer.colors))){
    pheno.pal <- cancer.colors[names(pheno.map)]
    pheno.names <- cancer.names[names(pheno.map)]
  }else{
    pheno.pal <- greyscale.palette(length(pheno.map)+2)[2:(length(pheno.map)+1)]
    pheno.names <- names(pheno.map)
    names(pheno.names) <- names(pheno.pal) <- names(pheno.map)
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)

  # Draw rectangles with labels
  y.start <- 0
  if(do.pop){
    rect(xleft=pop.starts, xright=pop.ends,
         ybottom=y.start + rect.hex.buffer,
         ytop=y.start + 1 - rect.hex.buffer,
         col=pop.pal[pop.labs], border=pop.pal[pop.labs])
    pop.rect.df <- data.frame("start" = pop.starts, "end" = pop.ends,
                              "len" = pop.ends-pop.starts, "lab" = pop.labs)
    pop.rect.df <- pop.rect.df[order(-pop.rect.df$len), ]
    sapply(unique(pop.labs), function(pop){
      pop.lab.x <- mean(as.numeric(pop.rect.df[which(pop.rect.df$lab == pop),
                                               c("start", "end")]))
      short.lab <- shorten.text(pop.names[pop],
                                width=pop.rect.df$len[which(pop.rect.df$lab == pop)],
                                cex=global.scaling.cex)
      text(x=pop.lab.x, y=y.start+0.5,
           labels=short.lab,
           col=optimize.label.color(pop.pal[pop], cutoff=0.8),
           cex=global.scaling.cex)
    })
    y.start <- y.start + 1
  }
  if(do.pheno){
    rect(xleft=pheno.starts, xright=pheno.ends,
         ybottom=y.start + rect.hex.buffer, ytop=y.start + 1 - rect.hex.buffer,
         col=pheno.pal[pheno.labs], border=pheno.pal[pheno.labs])
    pheno.rect.df <- data.frame("start" = pheno.starts, "end" = pheno.ends,
                                "len" = pheno.ends-pheno.starts, "lab" = pheno.labs)
    pheno.rect.df <- pheno.rect.df[order(-pheno.rect.df$len), ]
    sapply(unique(pheno.labs), function(pheno){
      pheno.lab.xs <- as.numeric(head(pheno.rect.df[which(pheno.rect.df$lab == pheno),
                                                        c("start", "end")], 1))
      text(x=mean(pheno.lab.xs), y=y.start+0.5,
           labels=shorten.text(pheno.names[pheno], abs(diff(pheno.lab.xs))),
           col=optimize.label.color(pheno.pal[pheno], cutoff=0.8),
           cex=global.scaling.cex)
    })
  }

  # Add labels
  sapply(1:length(group.labels), function(y){
    axis(2, at=y-0.5, tick=F, line=-1.2, labels=group.labels[y], las=2,
         cex.axis=global.scaling.cex)
  })
}

# Wrapper for main waterfall plot
main.waterfall <- function(gt.counts, out.prefix, pop=NULL, pheno=NULL,
                           pop.space.wex=0.025, pheno.space.wex=0.005){
  # Determine number of panels and figure sizing
  vcs <- intersect(names(var.class.names), unique(gt.counts$class))
  n.panels <- length(vcs)
  has.pop <- !is.null(pop)
  has.pheno <- !is.null(pheno)
  bottom.tracks <- sum(c(has.pop, has.pheno))
  has.bottom.tracks <- bottom.tracks > 0
  samples <- unique(gt.counts$sample)

  # Compute each sample's total variant counts
  total.counts <- as.data.frame(do.call("cbind", lapply(vcs, function(vc){
    sapply(samples, function(sid){
      s.v <- as.numeric(gt.counts[which(gt.counts$sample == sid
                                        & gt.counts$class == vc
                                        & is.na(gt.counts$subclass)), c("het", "hom")])
      if(length(s.v) == 0){0}else{sum(s.v, na.rm=T)}
    })
  })))
  colnames(total.counts) <- vcs
  rownames(total.counts) <- samples

  # Determine sample order as follows:
  # 1. Group samples on X axis by population from most to least total variants
  # 2. Order samples within each group first by phenotype (largest to smallest N)
  # 3. Order samples within `2` from most to least variation
  if(has.pop){
    pop <- pop[samples]
    pops <- unique(pop)
    pop.m <- sapply(pops, function(p){
      median(apply(total.counts[names(pop[which(pop == p)]), vcs, drop=FALSE], 1, sum))
    })
    pop.rank <- 1:length(pop.m)
    names(pop.rank) <- names(pop.m)[order(pop.m, decreasing=TRUE)]
  }
  if(has.pheno){
    pheno <- pheno[samples]
    pheno.k <- table(pheno)
    pheno.rank <- 1:length(pheno.k)
    names(pheno.rank) <- names(pheno.k)[order(pheno.k, decreasing=TRUE)]
  }
  s.z <- apply(apply(total.counts, 2, scale), 1, mean)
  s.rank <- 1:length(samples)
  names(s.rank) <- rownames(total.counts)[order(s.z, decreasing=T)]
  order.df <- data.frame("counts" = s.rank, row.names=names(s.rank))
  if(has.pop){
    order.df$pop <- pop.rank[pop[rownames(order.df)]]
  }
  if(has.pheno){
    order.df$pheno <- pheno.rank[pheno[rownames(order.df)]]
  }
  order.df <- order.df[, intersect(c("pop", "pheno", "counts"),
                                   colnames(order.df))]
  order.df <- order.df[do.call(order, order.df), ]
  samples.sorted <- rownames(order.df)

  # Define plot parameters
  # Target dimensions for 3 rows + 2 marker fields: 6.6" wide x 4" tall
  # Roughly 1.2" for each main row, 0.125" for each marker row, and 0.15" for bottom axis
  pdf.width <- 6.6
  panel.height <- 3.55/3
  marker.row.height <- 0.125
  bottom.axis.height <- 0.15
  pdf.height <- (panel.height*n.panels) + (bottom.tracks*marker.row.height) + bottom.axis.height
  layout.heights <- c(rep(panel.height, n.panels),
                      if(has.bottom.tracks){bottom.tracks*marker.row.height})
  layout.heights[length(layout.heights)] <- layout.heights[length(layout.heights)] + bottom.axis.height
  text.scaling.factor <- (7 / seq(10, 7, length.out=3))[n.panels]

  # Generate waterfall plot
  pdf(paste(out.prefix, "variants_per_genome.waterfall.pdf", sep="."),
      height=pdf.height, width=pdf.width)
  layout(matrix(1:(n.panels+has.bottom.tracks), byrow=T, ncol=1),
         heights=layout.heights)
  for(vc in vcs){
    waterfall.parmar <- c(0.25, 3.5, 0.5, 2.5) * text.scaling.factor
    if(!has.bottom.tracks & vc == vcs[n.panels]){
      parmar[1] <- 1.35 * text.scaling.factor
    }
    plot.count.waterfall(gt.counts, vc, pop, pheno, samples.sorted,
                         pop.space.wex, pheno.space.wex,
                         global.scaling.cex=text.scaling.factor,
                         parmar=waterfall.parmar)
  }
  if(has.bottom.tracks){
    add.waterfall.markers(order.df, pop.rank, pheno.rank,
                          pop.space.wex, pheno.space.wex, global.scaling.cex=text.scaling.factor,
                          parmar=c(1.2, 3.5, 0.1, 2.5)*text.scaling.factor)
  }
  mtext(paste(prettyNum(length(samples), big.mark=","), "germline genomes"),
        1, line=par("mar")[1]-1.1, cex=5/6)
  dev.off()
}

# Helper function to handle scatterplots used for inter-class comparisons
interclass.scatter <- function(plot.df, pop=NULL, title=NULL, label.units=NULL,
                               axis.lab.suffix=NULL, xlims=NULL, ylims=NULL,
                               diag.lm=F, cor.lm=F, stat.lab.hex=0.04,
                               pt.cex=0.65, parmar=c(2, 2.75, 1, 3.5)){
  # Get pointwise plotting parameters
  pw.params <- scatterplot.point.params(nrow(plot.df))
  if(is.null(pop)){
    pw.col <- rep("black", nrow(plot.df))
  }else{
    n.pops <- length(unique(pop))
    pop <- toupper(pop)
    if(all(pop %in% names(pop.colors))){
      pop.pal <- pop.colors[sort(unique(pop))]
    }else{
      pop.pal <- categorical.rainbow(n.pops)
      names() <- unique(pop)
    }
    pw.col <- pop.colors[pop[rownames(plot.df)]]
  }

  # Prepare plotting area
  if(is.null(xlims)){
    xlims <- range(plot.df[, 1])
  }
  if(is.null(ylims)){
    ylims <- range(plot.df[, 2])
  }
  prep.plot.area(xlims, ylims, parmar, xaxs="r", yaxs="r")
  if(diag.lm){abline(0, 1, col=annotation.color)}

  # Add points
  points(plot.df, pch=pw.params$pch, cex=0.75*pw.params$cex, xpd=T,
         col=sapply(pw.col, adjustcolor, alpha=pw.params$alpha))
  if(cor.lm){abline(lm(plot.df[, 2] ~ plot.df[, 1]), lty=2)}

  # Add axes & title
  mtext(3, text=title)
  ext.vc.abbrevs <- c("all" = "Total", var.class.abbrevs)
  clean.axis(1, label.units=label.units, infinite=T, min.ticks=2, max.ticks=4,
             title=tryCatch(paste(ext.vc.abbrevs[colnames(plot.df)[1]],
                                  axis.lab.suffix),
                            error=function(e){"X"}),
             label.line=-0.85, title.line=0, max.label.decimals=1)
  clean.axis(2, label.units=label.units, infinite=T, min.ticks=2, max.ticks=4,
             title=tryCatch(paste(ext.vc.abbrevs[colnames(plot.df)[2]],
                                  axis.lab.suffix),
                            error=function(e){"Y"}),
             label.line=-0.75, title.line=0.75, max.label.decimals=1)

  # Add title & legends
  r2 <- cor(plot.df[, 1], plot.df[, 2])^2
  ideal.quad <- head(setdiff(calc.plot.quadrant.density(plot.df[, 1], plot.df[, 2]),
                             c("I", "III")), 1)
  if(ideal.quad == "IV"){
    text(x=par("usr")[2]+(stat.lab.hex*diff(par("usr")[1:2])),
         y=par("usr")[3]+(2.5*stat.lab.hex*diff(par("usr")[3:4])),
         pos=2, cex=5/6, xpd=T,
         labels=bquote(italic(R)^2*"="*.(formatC(round(r2, 2), digits=2))))
  }else{
    text(x=par("usr")[1]-(stat.lab.hex*diff(par("usr")[1:2])),
         y=par("usr")[4]-(2.5*stat.lab.hex*diff(par("usr")[3:4])),
         pos=4, cex=5/6, xpd=T,
         labels=bquote(italic(R)^2*"="*.(formatC(round(r2, 2), digits=2))))
  }
  if(!is.null(pop)){
    pop.labs.at <- sort(sapply(names(pop.pal), function(p){
      mean(plot.df[names(pop)[which(pop == p)], 2], na.rm=T)
    }))
    if(all(pop %in% names(pop.abbreviations))){
      pop.labs <- pop.abbreviations[names(pop.labs.at)]
    }else{
      pop.labs <- names(pop.labs.at)
    }
    pop.labs.at <- smart.spacing(pop.labs.at, min.dist=0.025*diff(par("usr")[3:4]),
                                 lower.limit=par("usr")[3], upper.limit=par("usr")[4])
    names(pop.labs.at) <- names(pop.labs)
    yaxis.legend(pop.labs, x=par("usr")[2], y.positions=pop.labs.at,
                 sep.wex=0.05*diff(par("usr")[1:2]), colors=pop.pal[names(pop.labs)],
                 lwd=4, min.label.spacing=0.125*diff(par("usr")[3:4]))
  }

  return(r2)
}

# Scatterplot of sample heterozygosity between two classes of variants
plot.heterozygosity <- function(gt.counts, vc1, vc2, pop=NULL,
                                title="Heterozygosity"){
  # Compute heterozygosity for all samples for each variant class
  samples <- unique(gt.counts$sample)
  h.df <- as.data.frame(do.call("cbind", lapply(c(vc1, vc2), function(vc){
    sapply(samples, function(sid){
      sv1 <- gt.counts[(which(gt.counts$sample == sid
                              & gt.counts$class == vc
                              & is.na(gt.counts$subclass))),
                       c("het", "hom")]
      if(length(sv1) == 0){NA}else{sv1$het / sum(sv1)}
    })
  })))
  rownames(h.df) <- samples
  h.df <- h.df[which(complete.cases(h.df)), ]
  colnames(h.df) <- c(vc1, vc2)

  # Reorder samples from smallest to largest pop, if optioned
  if(!is.null(pop)){
    pop.k <- table(pop[rownames(h.df)])
    pop.idx <- length(pop.k):1
    names(pop.idx) <- names(pop.k)[order(pop.k)]
    h.df <- h.df[order(pop.idx[pop[rownames(h.df)]], h.df[, vc1], h.df[, vc2]), ]
  }

  # Generate plot
  xlims <- ylims <- range(h.df[, c(vc1, vc2)])
  r2 <- interclass.scatter(h.df, pop, title, label.units="percent",
                           axis.lab.suffix="het. rate", xlims=xlims,
                           ylims=ylims, diag.lm=T)

  return(c(r2, nrow(h.df)))
}

# Scatterplot of overall variant counts between two classes of variants
plot.count.comparisons <- function(gt.counts, vc1, vc2, pop=NULL,
                                   title="Variant count"){
  # Compute total number of variants for all samples for each variant class
  samples <- unique(gt.counts$sample)
  c.df <- as.data.frame(do.call("cbind", lapply(c(vc1, vc2), function(vc){
    sapply(samples, function(sid){
      sum(gt.counts[(which(gt.counts$sample == sid
                           & gt.counts$class == vc
                           & is.na(gt.counts$subclass))),
                    c("het", "hom")])
    })
  })))
  rownames(c.df) <- samples
  c.df <- c.df[which(complete.cases(c.df)), ]
  colnames(c.df) <- c(vc1, vc2)

  # Reorder samples from smallest to largest pop, if optioned
  if(!is.null(pop)){
    pop.k <- table(pop[rownames(c.df)])
    pop.idx <- length(pop.k):1
    names(pop.idx) <- names(pop.k)[order(pop.k)]
    c.df <- c.df[order(pop.idx[pop[rownames(c.df)]], c.df[, vc1], c.df[, vc2]), ]
  }

  # Generate plot
  r2 <- interclass.scatter(c.df, pop, title, label.units="count",
                           axis.lab.suffix="count", cor.lm=T)

  return(c(r2, nrow(c.df)))
}

# Wrapper function to handle inter-class comparison plots
plot.vc.comparisons <- function(gt.counts, vc1, vc2, out.prefix, pop=NULL){

  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.prefix <- paste(vc1, "vs", vc2, sep="_")

  # Scatterplot of heterozygosity per sample between all pairs of variant classes
  png(paste(out.prefix, ss.prefix, "heterozygosity.png", sep="."),
      height=2.25*300, width=2.75*300, res=300)
  m.tmp <- plot.heterozygosity(gt.counts, vc1, vc2, pop=pop)
  dev.off()

  ss.df[1, ] <- c(paste(ss.prefix, "heterozygosity_cor", sep="."), "r2", m.tmp)

  # Scatterplot of total variant counts per sample between all pairs of variant classes
  png(paste(out.prefix, ss.prefix, "variant_counts.png", sep="."),
      height=2.25*300, width=2.75*300, res=300)
  m.tmp <- plot.count.comparisons(gt.counts, vc1, vc2, pop=pop)
  dev.off()

  ss.df[2, ] <- c(paste(ss.prefix, "variant_count_cor", sep="."), "r2", m.tmp)

  return(ss.df)
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot variant distributions per genome")
parser$add_argument("--genotype-dist-tsv", required=TRUE, metavar=".tsv",
                    help=paste("Input .tsv of genotype counts per sample ",
                               "generated by clean_sample_genotypes.py"))
parser$add_argument("--ancestry-labels", metavar=".tsv",
                    help=paste("Optional two-column .tsv mapping sample IDs to ",
                               "ancestry labels. If provided, all plots will ",
                               "distinguish between ancestries."))
parser$add_argument("--phenotype-labels", metavar=".tsv",
                    help=paste("Optional two-column .tsv mapping sample IDs to ",
                               "phenotype labels. If provided, some plots will ",
                               "group samples by phenotype for easier comparisons."))
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
# args <- list("genotype_dist_tsv" = "~/scratch/dfci-g2c.v1.initial_qc.genotype_distribution.merged.tsv.gz",
#              "ancestry_labels" = "~/scratch/dfci-g2c.v1.qc_ancestry.tsv",
#              "phenotype_labels" = "~/scratch/dfci-g2c.v1.qc_phenotype.tsv",
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/g2c_vcf_qc_dev")

# # DEV (single variant class):
# args <- list("genotype_dist_tsv" = "~/scratch/dfci-ufc.sv.v1.initial_qc.genotype_distribution.merged.tsv.gz",
#              "ancestry_labels" = "~/scratch/dfci-g2c.v1.qc_ancestry.tsv",
#              "phenotype_labels" = "~/scratch/dfci-g2c.v1.qc_phenotype.tsv",
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/ufc_sv_qc_dev")

# # DEV (YL):
# args <- list("genotype_dist_tsv" = "~/scratch/YL.v1.2.genotype_distribution.merged.tsv.gz",
#              "ancestry_labels" = "~/scratch/YL_ancestry_labels.tsv",
#              "phenotype_labels" = "~/scratch/YL_pheno_labels.tsv",
#              "subset_samples" = "~/Desktop/Collins/VanAllen/jackie_younglung/younglung_metadata/YL.samples_with_complete_variant_data.list",
#              "custom_constants" = "~/Desktop/Collins/VanAllen/jackie_younglung/data/yl_qc.custom_constants.R",
#              "out_prefix" = "~/scratch/yl_terra_dev")

# Load custom constants if optioned
if(!is.null(args$custom_constants)){
  source(args$custom_constants)
}

# Load sample genotype counts
gt.counts <- load.gt.counts(args$genotype_dist_tsv)
samples <- unique(gt.counts$sample)

# Subset to samples in --subset-samples if optioned
if(!is.null(args$subset_samples)){
  keep.samples <- unique(read.table(args$subset_samples, header=F)[, 1])
  samples <- intersect(samples, keep.samples)
  gt.counts <- gt.counts[which(gt.counts$sample %in% keep.samples), ]
}

# Load sample ancestry and phenotypes, if optioned
pop <- load.labels(args$ancestry_labels, samples)
if(!is.null(pop)){
  gt.counts <- gt.counts[which(gt.counts$sample %in% names(pop)), ]
}
pheno <- load.labels(args$phenotype_labels, samples)
if(!is.null(pheno)){
  gt.counts <- gt.counts[which(gt.counts$sample %in% names(pheno)), ]
}

# Triple waterfall plot of counts per sample by class
main.waterfall(gt.counts, args$out_prefix, pop, pheno)

# Collect median counts per sample by class, subclass, frequency, zygosity, pop, and pheno
count.ss <- gather.count.sumstats(gt.counts, samples, pop, pheno)

# TODO: maybe add Sankey diagrams here

# Inter-class comparisions across samples
ic.ss <- lapply(list(c("all", "snv"), c("all", "indel"), c("all", "sv"),
                     c("snv", "indel"), c("snv", "sv"), c("indel", "sv")),
                function(vcs){
  vc1 <- vcs[1]; vc2 <- vcs[2]
  if(all(vcs %in% gt.counts$class)){
    plot.vc.comparisons(gt.counts, vc1, vc2, args$out_prefix, pop)
  }
})
if(length(ic.ss) > 0){
  ic.ss <- as.data.frame(do.call("rbind", ic.ss))
}
if(nrow(ic.ss) > 0){
  count.ss <- as.data.frame(rbind(count.ss, ic.ss))
}

# Write summary statistics to output file
colnames(count.ss)[1] <- paste("#", colnames(count.ss)[1], sep="")
write.table(count.ss, paste(args$out_prefix, "variants_per_sample.summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

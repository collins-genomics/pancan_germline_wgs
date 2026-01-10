#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Helper functions for VCF QC


#' Compute genotype benchmarking summary statistics
#'
#' Calculate summary metrics from output of [G2CR::get.gt.bench.plot.data()]
#'
#' @param plot.dat List of benchmarking data processed by [G2CR::get.gt.bench.plot.data()]
#' @param ss.prefix Character prefix for summary statistic descriptions
#' @param summarize.all Should a final series of summary rows be added to reflect
#' the overall summary statistics of all data in `plot.dat`? \[default: TRUE\]
#' @param trio.mode Is the input data from parent-child trios? \[default: FALSE\]
#'
#' @returns data.frame of summary statistics for various strata within `plot.dat`
#'
#' @param seealso [G2CR::get.gt.bench.plot.data()]
#'
#' @export calc.gt.bench.ss
#' @export
calc.gt.bench.ss <- function(plot.dat, ss.prefix, summarize.all=TRUE,
                             trio.mode=FALSE){
  ss.prefix <- sub("[ /]+", "_", ss.prefix)
  query.key <- if(trio.mode){"pass.rate"}else{"loose"}
  all.res <- lapply(1:length(plot.dat), function(i){
    prefix.i <- paste(ss.prefix, names(plot.dat)[i], sep=".")
    df.i <- plot.dat[[i]]

    # Gather marginal results within each evaluation set
    res.i <- as.data.frame(do.call("rbind", lapply(1:length(df.i), function(j){
      prefix.j <- paste(prefix.i, tolower(names(df.i)[j]), sep=".")
      df.j <- df.i[[j]]
      rbind(c(prefix.j, "median", df.j["Median", c(query.key, "n")]),
            c(prefix.j, "dynamic_range", dynamic.range(df.j[c("Min.", "Max."), query.key]), NA))
    })))
    colnames(res.i) <- c("analysis", "measure", "value", "n")

    # Get overall metric for each stratum across all evaluation sets
    i.med <- weighted.mean(as.numeric(sapply(df.i, function(d){d["Median", query.key]})),
                           as.numeric(sapply(df.i, function(d){d["Median", "n"]})))
    i.med.n <- sum(as.numeric(sapply(df.i, function(d){d["Median", "n"]})))
    i.min <- weighted.mean(as.numeric(sapply(df.i, function(d){d["Min.", query.key]})),
                           as.numeric(sapply(df.i, function(d){d["Min.", "n"]})))
    i.max <- weighted.mean(as.numeric(sapply(df.i, function(d){d["Max.", query.key]})),
                           as.numeric(sapply(df.i, function(d){d["Max.", "n"]})))
    i.dr <- i.max / i.min
    res.i.add <- as.data.frame(rbind(c(prefix.i, "median", i.med, i.med.n),
                                     c(prefix.i, "dynamic_range", i.dr, NA)))
    colnames(res.i.add) <- colnames(res.i)
    rbind(res.i, res.i.add)
  })

  # Calculate overall summary metrics for all variants, if optioned
  if(summarize.all){
    strata.n <- as.numeric(sapply(all.res, function(d){tail(d[which(d$measure == "median"), "n"], 1)}))
    all.med <- weighted.mean(as.numeric(sapply(all.res, function(d){tail(d[which(d$measure == "median"), "value"], 1)})),
                             strata.n)
    all.dr <- weighted.mean(as.numeric(sapply(all.res, function(d){tail(d[which(d$measure == "dynamic_range"), "value"], 1)})),
                            strata.n)
    res.all.add <- as.data.frame(rbind(c(ss.prefix, "median", all.med, sum(strata.n)),
                                       c(ss.prefix, "dynamic_range", all.dr, NA)))
    names(res.all.add) <- names(all.res[[1]])

    # Format and return sumstats as dataframe
    as.data.frame(rbind(do.call("rbind", all.res), res.all.add))
  }else{
    as.data.frame(do.call("rbind", all.res))
  }
}


#' Calculate genotype benchmarking statistics
#'
#' Calculate the benchmarking statistics necessary for plotting sample-level
#' genotype benchmarking
#'
#' @param bench.dat A list of data.frames of genotype benchmarking data loaded
#' with [G2CR::load.gt.benchmark.tsvs()]
#' @param vc Optional character vector of variant classes for subsetting
#' prior to calculations
#' @param vsc Optional character vector of variant subclasses for subsetting
#' prior to calculations
#' @param freq Optional character vector of frequency categories for subsetting
#' prior to calculations
#' @param zygosity Optional character vector of zygosities to consider for calculations
#' @param key.cols Numeric vector of column indexes to treat as
#' unique keys \[default: 1:5\]
#' @param trio.mode Is the input data from parent-child trios? \[default: FALSE\]
#'
#' @returns List of data.frames of benchmarking summary statistics amenable with
#' downstream plotting functions
#'
#' @seealso [G2CR::load.gt.benchmark.tsvs()]
#'
#' @export calc.gt.bench.stats
#' @export
calc.gt.bench.stats <- function(bench.dat, vc=NULL, vsc=NULL, freq=NULL,
                                zygosity=NULL, key.cols=NULL, trio.mode=FALSE){
  # Filter bench.dat according to vc, vsc, freq.breaks, and zygosity
  if(!is.null(vc)){
    bench.dat <- lapply(bench.dat,
                        function(bd.df){bd.df[which(bd.df$class %in% vc), ]})
  }
  if(!is.null(vsc)){
    bench.dat <- lapply(bench.dat,
                        function(bd.df){bd.df[which(bd.df$subclass %in% vsc), ]})
  }
  if(!is.null(freq)){
    bench.dat <- lapply(bench.dat,
                        function(bd.df){bd.df[which(bd.df$freq_bin %in% freq), ]})
  }
  if(!is.null(zygosity) & !trio.mode){
    bench.dat <- lapply(bench.dat,
                        function(bd.df){bd.df[which(bd.df$zygosity %in% zygosity), ]})
  }

  # Adjust key columns based on options
  if(is.null(key.cols)){
    key.cols <- if(trio.mode){1:4}else{1:5}
  }

  # Collapse each entry in bench.dat by sample
  sum.dat <- lapply(bench.dat, function(bd.df){
    s.key.name <- if(trio.mode){"family_id"}else{"sample"}
    samples <- as.character(unique(bd.df[, s.key.name]))
    s.df <- as.data.frame(do.call("rbind", lapply(samples, function(sid){
      apply(bd.df[which(bd.df[, s.key.name] == sid), -key.cols], 2, sum)
    })))
    rownames(s.df) <- samples
    s.df$n <- apply(s.df, 1, sum, na.rm=T)
    if(!trio.mode){
      if(nrow(s.df) > 0){
        s.df$strict <- s.df$gt_match / s.df$n
        s.df$loose <- apply(s.df[, c("gt_match", "carrier_match")],
                            1, sum, na.rm=T) / s.df$n
      }else{
        s.df <- data.frame("n" = numeric(),
                           "strict" = numeric(),
                           "loose" = numeric())
      }
      return(s.df[, c("n", "strict", "loose")])
    }else{
      if(nrow(s.df) > 0){
        s.df$pass.rate <- s.df$pass / s.df$n
      }else{
        s.df <- data.frame("n" = numeric(),
                           "pass.rate" =)
      }
      return(s.df[, c("n", "pass.rate")])
    }
  })

  # Summarize data for plotting
  lapply(sum.dat, function(sd.df){
    if(nrow(sd.df) > 0){
      apply(sd.df, 2, function(v){summary(v[which(!is.na(v) & !is.infinite(v))])})
    }else{
      summary.df <- as.data.frame(matrix(rep(NA, 6*ncol(sd.df)), nrow=6))
      colnames(summary.df) <- colnames(sd.df)
      rownames(summary.df) <- c("Min.", "1st Qu.", "Median", "Mean",
                                "3rd Qu.", "Max.")
      return(summary.df)
    }
  })
}


#' Calculate variants per genome for genotype benchmarking
#'
#' Compute total variants per genome from all elements in the output of [G2CR::get.gt.bench.plot.data()]
#'
#' @param plot.dat List of benchmarking data processed by [G2CR::get.gt.bench.plot.data()]
#'
#' @returns integer of the average number of variants per genome
#'
#' @export calc.gt.bench.vpg
#' @export
calc.gt.bench.vpg <- function(plot.dat){
  sum(unlist(sapply(plot.dat, function(lo){
    sapply(lo, function(li){li["Median", "n"]})
  })), na.rm=T)
}


#' Clean size or AF break labels
#'
#' Converts size or AF breaks encoded from characters to numeric values
#'
#' @param vals A character vector of break values to be parsed
#'
#' @returns A named vector of numeric values corresponding to `vals`, where
#' element names correspond to the original entries in `vals`
#'
#' @export clean.breaks
#' @export
clean.breaks <- function(vals){
  orig <- vals
  vals <- gsub("[eE]-", " * 10 ^ -", vals)
  vals <- gsub("[eE]\\+", " * 10 ^ ", vals)
  vals <- gsub("[A-Z]|[a-z]|_", "", vals)
  vals <- sapply(vals, function(v){as.numeric(eval(parse(text=v)))})
  names(vals) <- orig
  return(vals)
}


#' Gather genotype benchmarking plot data
#'
#' Wrapper to handle collection & formatting of data for plotting genotype benchmarking
#'
#' @param bench.dat List of benchmarking data.frames processed by [G2CR::load.gt.benchmark.tsvs()]
#' @param vc Character vector of variant classes for subsetting. See [G2CR::calc.gt.bench.stats()]
#' @param trio.mode Is the input data from parent-child trios? \[default: FALSE\]
#'
#' @returns List of data.frames, one each for rare heterozygous genotypes, common
#' heterozygous genotypes, and all homozygous genotypes
#'
#' @seealso [G2CR::load.gt.benchmark.tsvs()], [G2CR::calc.gt.bench.stats()]
#'
#' @export get.gt.bench.plot.data
#' @export
get.gt.bench.plot.data <- function(bench.dat, vc=NULL, trio.mode=FALSE){
  # Check frequency bins in the data and order from rare -> common
  freq.breaks <- clean.breaks(unique(bench.dat[[1]]$freq_bin))
  freq.breaks <- freq.breaks[order(-freq.breaks, names(freq.breaks),
                                   decreasing=T)]

  # By default, we want to split by rare het, common het, and all homozygotes
  # We ignore zygosity for trios as that is less simple to define across three genomes
  if(!trio.mode){
    list("rare.het"=calc.gt.bench.stats(bench.dat, vc=vc,
                                        freq=head(names(freq.breaks), 1),
                                        zygosity="het"),
         "common.het"=calc.gt.bench.stats(bench.dat, vc=vc,
                                          freq=tail(names(freq.breaks), 1),
                                          zygosity="het"),
         "hom"=calc.gt.bench.stats(bench.dat, vc=vc, zygosity="hom"))
  }else{
    list("rare"=calc.gt.bench.stats(bench.dat, vc=vc,
                                    freq=head(names(freq.breaks), 1),
                                    trio.mode=TRUE),
         "common"=calc.gt.bench.stats(bench.dat, vc=vc,
                                      freq=tail(names(freq.breaks), 1),
                                      trio.mode=TRUE))
  }
}


#' Load genotype benchmarking data
#'
#' Load & postprocess a set of genotype benchmarking .tsvs
#'
#' @param tsvs.in Character vector of paths to .tsv inputs
#' @param set.names Character names for each .tsv in `tsvs.in`
#' @param key.cols Numeric vector of column indexes to treat as
#' unique keys \[default: 1:5\]
#' @param trio.mode Is the input data from parent-child trios? \[default: FALSE\]
#' @param keep.samples Optional list of sample IDs to retain from `tsvs.in`. By
#' default, all samples will be retained
#' @param keep.vcs Optional list of variant classes to retain from `tsvs.in`. By
#' default, all variant classes will be retained
#'
#' @returns List of one data.frame per element in `tsvs.in`
#'
#' @export load.gt.benchmark.tsvs
#' @export
load.gt.benchmark.tsvs <- function(tsvs.in, set.names, key.cols=1:5,
                                   trio.mode=FALSE, keep.samples=NULL,
                                   keep.vcs=NULL){
  # Do nothing if no .tsv files are provided
  if(is.null(tsvs.in) | length(tsvs.in) == 0){
    return(NULL)
  }

  # Error out if set.names doesn't match tsvs.in
  if(length(tsvs.in) != length(set.names)){
    stop("mismatch between --set-name and input .tsvs; check your inputs")
  }

  # Otherwise, load each .tsv as raw data.frame
  bench.dat <- lapply(tsvs.in, function(tsv.in){
    d <- read.table(tsv.in, sep="\t", comment.char="", check.names=F, header=T)
    colnames(d)[1] <- gsub("#", "", colnames(d)[1], fixed=T)
    if(!is.null(keep.samples)){
      d <- d[which(d[, 1] %in% keep.samples), ]
    }
    if(!is.null(keep.vcs)){
      d <- d[which(d[, 2] %in% keep.vcs), ]
    }
    return(d)
  })
  names(bench.dat) <- set.names

  # Get union of all samples and row metadata
  s.key <- if(trio.mode){"family_id"}else{"sample"}
  samples <- unique(unlist(sapply(bench.dat, function(d){d[, s.key]})))
  bench.strata <- unique(do.call("rbind", bench.dat)[, setdiff(key.cols, 1)])

  # Fill zeros for missing samples in any set or strata to ensure consistent
  # dimensionality across all strata and datasets
  lapply(bench.dat, function(bd.df){
    add.df <- do.call("rbind", lapply(1:nrow(bench.strata), function(si){
      s.hits <- merge(bd.df, bench.strata[si, ], all.y=T, sort=F)[, s.key]
      s.miss <- setdiff(samples, s.hits)
      if(length(s.miss) == 0){
        return(NULL)
      }else{
        add.df <- data.frame(s.miss)
        colnames(add.df) <- s.key
        for(sx in 1:ncol(bench.strata)){
          add.df[, colnames(bench.strata)[sx]] <- bench.strata[si, sx]
        }
        for(oi in setdiff(1:ncol(bd.df), key.cols)){
          add.df[, colnames(bd.df)[oi]] <- 0
        }
        return(add.df)
      }
    }))
    bd.df <- as.data.frame(rbind(bd.df, add.df))

    # For sake of plotting, collapse complex SVs, inversions, and translocations
    # This is a practical consideration because simple inversions and translocations
    # are super sparse and therefore don't need to be treated separately
    if(length(intersect(c("INV", "CTX"), bd.df$subclass)) > 0){
      cpx.rows <- which(bd.df$subclass %in% c("CPX", "INV", "CTX"))
      cpx.df <- bd.df[cpx.rows, ]
      nocpx.df <- bd.df[-cpx.rows, ]
      cpx.strata <- as.data.frame(unique(bench.strata[which(bench.strata$subclass %in% c("CPX", "INV", "CTX")), -(1:2)]))
      colnames(cpx.strata) <- colnames(bench.strata[-c(1:2)])
      cpx.add <- as.data.frame(do.call("rbind", lapply(1:nrow(cpx.strata), function(si){
        si.df <- merge(cpx.df, cpx.strata[si, ], sort=F)[, colnames(bd.df)]
        do.call("rbind", lapply(unique(si.df[, s.key]), function(sid){
          s.sub <- si.df[which(si.df[, s.key] == sid), -key.cols]
          if(length(s.sub) > 0){
            c(sid, "sv", "CPX", as.character(cpx.strata[si, ]), apply(s.sub, 2, sum))
          }else{
            c(sid, "sv", "CPX", as.character(cpx.strata[si, ]), rep(0, 3))
          }
        }))
      })))
      colnames(cpx.add) <- colnames(nocpx.df)
      bd.df <- as.data.frame(rbind(nocpx.df, cpx.add))
    }
    bd.df[, -key.cols] <- apply(bd.df[, -key.cols], 2, as.numeric)
    return(bd.df)
  })
}


#' Plot all genotype benchmarking comparisons
#'
#' Main wrapper to plot a single comparison (PPV or sensitivity) across all strata
#'
#' @param bench.dat A list of data.frames of genotype benchmarking data loaded
#' with [G2CR::load.gt.benchmark.tsvs()]
#' @param out.prefix Path and filename prefix for all plots
#' @param set.colors Vector of colors for each sub-element in `plot.dat`
#' \[default: greyscale\]
#' @param ref.title Character name for reference/external benchmarking dataset
#' \[default: "external benchmarking"\]
#' @param metric.name Name of metric being evaluated
#' @param title.preps Prepositions to add in title \[default: c("for", "in")\]
#'
#' @seealso [G2CR::load.gt.benchmark.tsvs()], [G2CR::plot.gt.bench()]
#'
#' @export plot.all.gt.bench.strata
#' @export
plot.all.gt.bench.strata <- function(bench.dat, out.prefix, set.colors,
                                     ref.title="external benchmarking", metric.name=NULL,
                                     title.preps=c("for", "in")){
  # Check what vcs are present in data
  vcs <- unique(as.character(sapply(bench.dat, function(d){d$class})))

  # Get unique number of samples across all strata
  n.samples <- length(unique(unlist(lapply(bench.dat, function(d){d$sample}))))

  # Get AF cutoff for common vs. rare
  af.cutoff <- min(unique(clean.breaks(unique(bench.dat[[1]]$freq_bin))))

  # Prepare collection dataframe for summary stat logging
  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.prefix <- sub("[ /]+", "_", paste("gt_benchmarking", tolower(ref.title), sep="."))

  # Plot all variants
  plot.dat <- get.gt.bench.plot.data(bench.dat)
  ss.df <- rbind(ss.df, calc.gt.bench.ss(plot.dat, paste(ss.prefix, "all", sep=".")))
  vpg <- calc.gt.bench.vpg(plot.dat)
  plot.title <- paste(metric.name, title.preps[1], clean.numeric.labels(vpg),
                      "variants/genome\n", title.preps[2],
                      clean.numeric.labels(n.samples), ref.title, "samples")
  pdf(paste(out.prefix, "all.barplot.pdf", sep="."), height=2.25, width=2.45)
  plot.gt.bench(plot.dat=plot.dat,
                strata.names=c("Het.", "Het.", "Hom."),
                set.colors=set.colors,
                metric.name=metric.name,
                title=plot.title,
                af.cutoff=af.cutoff)
  dev.off()

  # Plot each variant class individually, if present
  for(vc in vcs){
    plot.dat <- get.gt.bench.plot.data(bench.dat, vc=vc)
    ss.df <- rbind(ss.df, calc.gt.bench.ss(plot.dat, paste(ss.prefix, vc, sep=".")))
    vpg <- calc.gt.bench.vpg(plot.dat)
    plot.title <- paste(metric.name, title.preps[1],
                        clean.numeric.labels(vpg, min.label.length=2, acceptable.decimals=1),
                        gsub("I", "i", paste(var.class.abbrevs[vc], "s", sep="")),
                        "/ genome\n", title.preps[2],
                        clean.numeric.labels(n.samples), ref.title, "samples")
    if(nchar(metric.name) >= 10){
      plot.title <- gsub(" / ", "/", plot.title, fixed=T)
    }
    pdf(paste(out.prefix, vc, "barplot.pdf", sep="."), height=2.25, width=2.45)
    plot.gt.bench(plot.dat=plot.dat,
                  strata.names=c("Het.", "Het.", "Hom."),
                  set.colors=set.colors,
                  metric.name=metric.name,
                  title=plot.title,
                  af.cutoff=af.cutoff)
    dev.off()
  }

  # Plot short variants by subclass, if present
  if(length(intersect(c("snv", "indel"), vcs)) > 1){
    short.vscs <- unique(as.character(sapply(bench.dat, function(d){
      d$subclass[which(d$class %in% c("snv", "indel"))]
    })))

    plot.dat <- lapply(short.vscs, function(vsc){
      calc.gt.bench.stats(bench.dat, vsc=vsc)
    })
    names(plot.dat) <- short.vscs
    short.ss.df <- calc.gt.bench.ss(plot.dat, ss.prefix, summarize.all=FALSE)
    ss.df <- rbind(ss.df, short.ss.df)
    vpg <- calc.gt.bench.vpg(plot.dat)
    plot.title <- paste(metric.name, title.preps[1],
                        clean.numeric.labels(vpg, min.label.length=2, acceptable.decimals=1),
                        "short variants per\ngenome", title.preps[2],
                        clean.numeric.labels(n.samples), ref.title, "samples")
    pdf(paste(out.prefix, "short_variants.by_subclass.barplot.pdf", sep="."),
        height=2.25, width=1.05+(1.3*length(short.vscs)/3))
    plot.gt.bench(plot.dat=plot.dat,
                  strata.names=var.subclass.abbrevs[short.vscs],
                  set.colors=set.colors,
                  metric.name=metric.name,
                  title=plot.title)
    dev.off()
  }

  # Plot SVs by subclass, if present
  if("sv" %in% vcs){
    sv.vscs <- unique(as.character(sapply(bench.dat, function(d){
      setdiff(d$subclass[which(d$class =="sv")], "CNV")
    })))

    plot.dat <- lapply(sv.vscs, function(vsc){
      calc.gt.bench.stats(bench.dat, vsc=vsc)
    })
    names(plot.dat) <- sv.vscs
    sv.ss.df <- calc.gt.bench.ss(plot.dat, ss.prefix, summarize.all=FALSE)
    ss.df <- rbind(ss.df, sv.ss.df)
    vpg <- calc.gt.bench.vpg(plot.dat)
    plot.title <- paste(metric.name, title.preps[1],
                        clean.numeric.labels(vpg, min.label.length=2, acceptable.decimals=1),
                        "SVs/genome\n", title.preps[2],
                        clean.numeric.labels(n.samples), ref.title, "samples")
    pdf(paste(out.prefix, "SVs.by_subclass.barplot.pdf", sep="."),
        height=2.25, width=1.05+(1.4*length(sv.vscs)/3))
    plot.gt.bench(plot.dat=plot.dat,
                  strata.names=var.subclass.abbrevs[sv.vscs],
                  set.colors=set.colors,
                  metric.name=metric.name,
                  title=plot.title)
    dev.off()
  }

  # Return summary statistics for logging
  med.ridx <- which(ss.df$measure == "median")
  ss.df$measure[med.ridx] <- sub("[ /]+", "_",
                                 paste(ss.df$measure[med.ridx],
                                       tolower(metric.name), sep="_"))
  ss.df$measure[-med.ridx] <- sub("[ /]+", "_",
                                  paste(tolower(metric.name),
                                        ss.df$measure[-med.ridx], sep="_"))
  return(ss.df)
}


#' Plot genotype benchmarking results
#'
#' Plot a user-specified set of sample-level genotype benchmarking strata
#'
#' @param plot.dat List of genotype benchmarking data.frames loaded by
#' [G2CR::get.gt.bench.plot.data()]
#' @param strata.names Vector of names for each element in `plot.dat`
#' \[default: use `names(plot.dat)`\]
#' @param set.colors Vector of colors for each sub-element in `plot.dat`
#' \[default: greyscale\]
#' @param metric.name Name of metric being evaluated
#' @param title Title for plot
#' @param af.cutoff Allele frequency cutoff to distinguish between common and
#' rare variants; will be used to add second row to bottom X-axis with AF bins
#' @param strata.space Space between strata specified in user X units \[default: 0.15\]
#' @param bar.hex Relative horizontal expansion parameter for the median bars
#' and the outer "whiskers"; must be bounded on \[0, 1\] \[default: 0.6\]
#' @param whisker.hex Relative horizontal expansion for median whiskers \[default: 0.035\]
#' @param trio.mode Is the input data from parent-child trios? \[default: FALSE\]
#' @param parmar Value of `mar` passed to `par()`
#'
#' @export plot.gt.bench
#' @export
plot.gt.bench <- function(plot.dat, strata.names=NULL, set.colors=NULL,
                          metric.name=NULL, title=NULL, af.cutoff=NULL,
                          strata.space=0.15, bar.hex=0.6,
                          whisker.hex=0.035, trio.mode=FALSE,
                          parmar=c(1.8, 2.55, 1.8, 2.5)){
  # Drop strata that are strictly all NAs / have no data
  keep.strata <- which(sapply(plot.dat, function(pd.l){
    !all(is.na(do.call("rbind", pd.l)))
    }))
  if(length(keep.strata) == 0){
    return()
  }
  plot.dat <- plot.dat[keep.strata]

  # Get plot parameters
  n.strata <- length(plot.dat)
  n.sets <- length(plot.dat[[1]])
  sets.left <- which(round((1:n.sets) / n.sets) == 0)
  sets.right <- which(round((1:n.sets) / n.sets) == 1)
  if(is.null(strata.names)){
    strata.names <- names(plot.dat)
  }else{
    strata.names <- strata.names[keep.strata]
  }
  if(is.null(set.colors)){
    set.colors <- greyscale.palette(n.strata)
  }
  if(trio.mode){
    set.pals <- set.colors
  }else{
    set.pals <- lapply(set.colors, function(col.base){
      colorRampPalette(c("black", col.base, "white"))(5)[-c(1:2, 5)]
    })
  }
  whisker.key <- if(trio.mode){"pass.rate"}else{"loose"}
  iqr.legend.y <- if(trio.mode){c(0.2, 0.35, 0.5, 0.35)}else{c(0, 0.1, 0.2, 0.1)}

  # Get generic spacing for elements in each stratum
  s.x.start <- strata.space / 2
  s.x.end <- 1 - (strata.space / 2)
  s.width <- s.x.end - s.x.start
  s.bar.width <- s.width * bar.hex
  s.whisker.width <- s.width * (1 - bar.hex)
  whiskers.left <- length(sets.left)
  whiskers.right <- length(sets.right)
  whiskers.left.width <- s.whisker.width * whiskers.left / n.sets
  whiskers.right.width <- s.whisker.width * whiskers.right / n.sets
  s.bar.x.left <- s.x.start + whiskers.left.width
  s.bar.x.right <- s.bar.x.left + s.bar.width
  whiskers.left.x.at <- if(whiskers.left == 0){c()}else{
    seq(s.x.start, s.bar.x.left,
        length.out=whiskers.left+2)[-c(1, whiskers.left+2)]
  }
  whiskers.right.x.at <- if(whiskers.right == 0){c()}else{
    seq(s.bar.x.right, s.x.end,
        length.out=whiskers.right+2)[-c(1, whiskers.right+2)]
  }
  whisker.x.adj <- whisker.hex * diff(par("usr")[1:2])

  # Prep plot area
  prep.plot.area(xlims=c(0, n.strata), ylims=c(0, 1.05), parmar=parmar)
  segments(x0=par("usr")[1], x1=par("usr")[2], y0=1, y1=1,
           col=annotation.color, xpd=T, lty=3)

  # Add strata labels
  sapply(1:n.strata, function(x){
    axis(1, at=x-0.5, tick=F, line=-1, labels=strata.names[x],
         cex.axis=5.5/6)
  })
  if(!is.null(af.cutoff)){
    mtext("AF", side=1, at=par("usr")[1]-(0.0125*diff(par("usr")[1:2])),
          adj=1, line=0.8, cex=5/6, col=annotation.color)
    axis(1, at=0.5, tick=F, line=-0.2, cex.axis=5/6, col.axis=annotation.color,
         labels=bquote("" <= .(paste(100 * af.cutoff, "% ", sep=""))))
    axis(1, at=1.5, tick=F, line=-0.2, cex.axis=5/6, col.axis=annotation.color,
         labels=bquote("" > .(paste(100 * af.cutoff, "% ", sep=""))))
    axis(1, at=2.5, tick=F, line=-0.2, cex.axis=5/6, col.axis=annotation.color,
         labels="All")
  }

  # Add left Y axis
  clean.axis(2, at=seq(0, 1, 0.25), label.units="percent", max.ticks=6,
             title=metric.name, title.line=0.7)

  # Add upper title
  mtext(3, text=title, line=0, cex=5.5/6)

  # Add legend to right Y axis
  if(!trio.mode){
    right.y.cs <- c(0, as.numeric(plot.dat[[n.strata]][[n.sets]]["Median", c("strict", "loose")]))
    right.y.mids <- sapply(1:n.sets, function(k){mean(right.y.cs[k+(0:1)])})
    right.y.lab.pos <- yaxis.legend(c("GT", "Carrier"),
                                    x=s.bar.x.right+n.strata-1,
                                    y.positions=right.y.mids,
                                    sep.wex=0.05*diff(par("usr")[1:2]),
                                    upper.limit=0.55, lower.limit=0.35,
                                    colors=hex2grey(set.pals[[n.sets]]),
                                    label.cex=4.5/6, lwd=2, return.label.pos=T)
    axis(4, at=max(right.y.lab.pos)+0.1125, las=2, tick=F, line=-0.7,
         cex.axis=5/6, labels="Match:")
  }
  set.labs.y <- seq(0.95, 0.7, length.out=n.sets+1)[1:n.sets]
  points(x=rep(par("usr")[2]+strata.space, 2), y=set.labs.y,
         pch=15, col=set.colors, xpd=T)
  text(x=rep(par("usr")[2]+(strata.space/2), 2), y=set.labs.y-0.01,
       cex=5/6, xpd=T, labels=names(plot.dat[[1]]), pos=4)
  polygon(x=par("usr")[2]+strata.space+c(0, -0.03*diff(par("usr")[1:2]),
                                         0, 0.03*diff(par("usr")[1:2])),
          y=iqr.legend.y, col=adjustcolor(annotation.color, alpha=0.4),
          border=NA, xpd=T)
  points(x=rep(par("usr")[2]+strata.space, 3), y=iqr.legend.y[1:3],
         pch=c(25, 3, 24), bg=c(annotation.color, NA, annotation.color),
         col=c(NA, annotation.color, NA), lwd=c(0, 2, 0), xpd=T, cex=0.8)
  text(x=rep(par("usr")[2]+(strata.space/2), 3), y=iqr.legend.y[1:3],
       cex=4.5/6, xpd=T, labels=c("Min", "IQR", "Max"), pos=4)

  # Plot each stratum
  sapply(1:n.strata, function(i){
    x.inc <- i-1
    s.dat <- plot.dat[[i]]

    if(all(is.na(s.dat))){
      return(NULL)
    }

    # Draw left whiskers
    sapply(1:whiskers.left, function(k){
      k.x <- whiskers.left.x.at[k] + x.inc
      k.y <- as.numeric(s.dat[[sets.left[k]]][c("Min.", "1st Qu.", "Median",
                                                "3rd Qu.", "Max."), whisker.key])
      for(color in c("white", adjustcolor(set.colors[sets.left[k]], alpha=1/3))){
        polygon(x=c(k.x, k.x-whisker.x.adj, k.x, k.x+whisker.x.adj),
                y=k.y[c(1, 3, 5, 3)], border=NA, col=color)
      }
      segments(x0=k.x, x1=k.x, y0=k.y[2], y1=k.y[4],
               col=set.colors[sets.left[k]])
      segments(x0=k.x-whisker.x.adj, x1=k.x+whisker.x.adj,
               y0=k.y[3], y1=k.y[3], col=set.colors[sets.left[k]])
      points(x=c(k.x, k.x), y=k.y[c(1, 5)], bg=set.colors[sets.left[k]],
             pch=c(25, 24), xpd=T, col=NA, cex=0.4)
    })

    # Draw right whiskers
    sapply(1:whiskers.right, function(k){
      k.x <- whiskers.right.x.at[k] + x.inc
      k.y <- as.numeric(s.dat[[sets.right[k]]][c("Min.", "1st Qu.", "Median",
                                                 "3rd Qu.", "Max."), whisker.key])
      for(color in c("white", adjustcolor(set.colors[sets.right[k]], alpha=1/3))){
        polygon(x=c(k.x, k.x-whisker.x.adj, k.x, k.x+whisker.x.adj),
                y=k.y[c(1, 3, 5, 3)], border=NA, col=color)
      }
      segments(x0=k.x, x1=k.x, y0=k.y[2], y1=k.y[4], col=set.colors[sets.right[k]])
      segments(x0=k.x-whisker.x.adj, x1=k.x+whisker.x.adj,
               y0=k.y[3], y1=k.y[3], col=set.colors[sets.right[k]])
      points(x=c(k.x, k.x), y=k.y[c(1, 5)], bg=set.colors[sets.right[k]],
             pch=c(25, 24), xpd=T, col=NA, cex=0.4)
    })

    # Draw rectangles
    rect.w.n <- as.numeric(sapply(s.dat, function(d){d["Median", "n"]})) + 0.1
    rect.x <- quantile(c(s.bar.x.left, s.bar.x.right) + x.inc,
                       probs=c(0, cumsum(rect.w.n / sum(rect.w.n))))
    sapply(1:length(s.dat), function(k){
      if(trio.mode){
        ybottoms <- 0
        ytops <- as.numeric(s.dat[[k]]["Median", "pass.rate"])
      }else{
        ybottoms <- as.numeric(c(0, s.dat[[k]]["Median", "strict"]))
        ytops <- as.numeric(s.dat[[k]]["Median", c("strict", "loose")])
      }
      rect(xleft=rect.x[k], xright=rect.x[k+1], ybottom=ybottoms, ytop=ytops,
           border=set.pals[[k]], col=set.pals[[k]], lwd=0.5)
    })
  })
}


#' Import a compressed summary distribution
#'
#' Import a precomputed compressed summary distribution from a .tsv
#'
#' @param tsv.in Path to input .tsv with compressed distribution
#' @param key.cols Column indexes to be treated as categorical key values \[default: 1-2\]
#'
#' @returns A two-element list of `df`, a dataframe of the distribution encoded
#' in `tsv.in`, and `breaks`, a numeric vector corresponding to the column-axis
#' break values
#'
#' @seealso [G2CR::clean.breaks]
#'
#' @export read.compressed.distrib
#' @export
read.compressed.distrib <- function(tsv.in, key.cols=1:2){
  if(is.null(tsv.in)){return(NULL)}
  df <- read.table(tsv.in, header=T, sep="\t", comment.char="",
                   quote="", check.names=F)
  colnames(df)[1] <- gsub("#", "", colnames(df)[1])
  df[, -key.cols] <- apply(df[, -key.cols], 2, as.numeric)
  dist.breaks <- clean.breaks(colnames(df)[-key.cols])
  list("df" = df, "breaks" = dist.breaks)
}


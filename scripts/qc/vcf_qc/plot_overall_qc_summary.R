#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot overall QC summary figures


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(DescTools, quietly=TRUE)
require(G2CR, quietly=TRUE)
load.constants("all")

# Declare global constants
# List of metrics to transform to percentages when plotting
pct.metrics <- c("hwe_pass", "hwe_acc", "ld", "site_ratio", "site_sens",
                 "site_ppv", "trio_inh_rate", "rep_match_rate", "heterozygosity")
# List of metrics with undefined default targets
no.target.default <- c("site_count", "site_count.rare", "site_count.singletons",
                       "variants_per_genome")
target.map <- c("site_ratios.all:NA" = NA,
                "site_ratios.snv:ti_tv_ratio" = 2.1 / 3.1,
                "site_ratios.indel:ins_del_ratio" = 1 / 2,
                "site_ratios.sv:DUP_DEL_ratio" = 1 / 2)
for(vc in c("all", names(var.class.abbrevs))){
  # Set VC-uniform
  target.map[paste(vc, "common_hwe:pct_pass", sep=".")] <- 0.99
  target.map[paste(vc, "common_hwe:weighted_accuracy", sep=".")] <- 0.99
  target.map[paste("heterozygosity.", vc, ":median", sep="")] <- 1.55 / 2.55
  target.map[paste("heterozygosity.", vc, ":reciprocal_dynamic_range", sep="")] <- 0.8
  target.map[paste("variants_per_genome.", vc, ":reciprocal_dynamic_range", sep="")] <- 0.85

  # Set undefined
  target.map[paste(paste(setdiff(no.target.default, "variants_per_genome"),
                         vc, sep="."), "count", sep=":")] <- NA
  target.map[paste(paste("variants_per_genome", vc, sep="."), "median", sep=":")] <- NA
}


##################
# Data Functions #
##################
# Invert dynamic range entries
invert.dynamic.ranges <- function(ss){
  dr.idxs <- grep("dynamic_range", ss$measure)
  ss$value[dr.idxs] <- 1 / ss$value[dr.idxs]
  ss$measure[dr.idxs] <- gsub("dynamic_range", "reciprocal_dynamic_range", ss$measure[dr.idxs])
  return(ss)
}

# Format count stats for a single VC
get.counts <- function(ss, vc){
  # Gather data
  count.suffix <- vc
  prop.rare <- ss[which(ss$analysis == paste("pct_rare", vc, sep=".")), ]
  rare.k <- if(nrow(prop.rare) > 0){round(prop.rare$value * prop.rare$n, 0)}else{0}
  rare.n <- if(nrow(prop.rare) > 0){prop.rare$n}else{0}
  prop.stn <- ss[which(ss$analysis == paste("pct_singletons", vc, sep=".")), ]
  stn.k <- if(nrow(prop.stn) > 0){round(prop.stn$value * prop.stn$n, 0)}else{0}
  stn.n <- if(nrow(prop.stn) > 0){prop.stn$n}else{0}
  count.df <- as.data.frame(rbind(
    ss[which(ss$analysis == paste("site_count", count.suffix, sep=".")), ],
    c(paste("site_count.rare", count.suffix, sep="."), "count", rare.k, rare.n),
    c(paste("site_count.singletons", count.suffix, sep="."), "count", stn.k, stn.n),
    ss[which(ss$analysis == paste("variants_per_genome", vc, sep=".")
             & ss$measure == "median"), ]
  ))

  # Assign universal row names
  rownames(count.df) <-
    unlist(sapply(apply(count.df[, c("analysis", "measure")], 1, paste, collapse="."),
                  function(qstr){
                    if(qstr == paste(paste("site_count", count.suffix, "count", sep="."))){
                      "count_all"
                    }else if(startsWith(qstr, "site_count.rare")){
                      "count_rare"
                    }else if(startsWith(qstr, "site_count.singletons")){
                      "count_singleton"
                    }else if(startsWith(qstr, "variants_per_genome")){
                      "count_per_genome"
                    }
                  }))

  return(count.df)
}

# Format site benchmarking stats for a single VC
get.sb <- function(ss, vc, ref.prefix=NULL){
  # Collect data
  sb.df <- as.data.frame(rbind(
    ss[which(ss$analysis == paste(vc, "common_hwe", sep=".")
             & ss$measure == "pct_pass"), ],
    ss[which(ss$analysis == paste(vc, "common_hwe", sep=".")
             & ss$measure == "weighted_accuracy"), ],
    ss[which(ss$analysis == paste(vc, "common_ld.any", sep=".")
             & ss$measure == "tag_rate"), ],
    if(vc == "all"){
      data.frame("analysis" = "site_ratios.all", "measure" = NA, "value" = NA, "n" = NA)
    }else{
      ss[which(ss$analysis == paste("site_ratios", vc, sep=".")), ]
    },
    ss[which(ss$analysis == paste(ref.prefix, vc, "common_af_cor", sep=".")), ],
    ss[which(ss$analysis == gsub("^\\.|\\.$", "",
                                 paste(ref.prefix, "sensitivity.common",
                                       if(vc!="all"){vc}, sep="."))), ],
    ss[which(ss$analysis == gsub("^\\.|\\.$", "",
                                 paste(ref.prefix, "ppv.common",
                                       if(vc!="all"){vc}, sep="."))), ]
  ))

  # Assign universal row names
  rownames(sb.df) <-
    unlist(sapply(apply(sb.df[, c("analysis", "measure")], 1, paste, collapse="."),
                  function(qstr){
                    if(qstr == paste(paste(vc, "common_hwe.pct_pass", sep="."))){
                      "hwe_pass"
                    }else if(qstr == paste(paste(vc, "common_hwe.weighted_accuracy", sep="."))){
                      "hwe_acc"
                    }else if(qstr == paste(paste(vc, "common_ld.any.tag_rate", sep="."))){
                      "ld"
                    }else if(startsWith(qstr, "site_ratios")){
                      "site_ratio"
                    }else if(endsWith(qstr, "common_af_cor.r2")){
                      "site_af_cor"
                    }else if(endsWith(qstr, ".sensitivity")){
                      "site_sens"
                    }else if(endsWith(qstr, ".ppv")){
                      "site_ppv"
                    }
                  }))

  return(sb.df)
}

# Format genotype benchmarking stats for a single VC
get.gb <- function(ss, vc, gb.prefixes=c()){
  # Collect data
  gb.df <- as.data.frame(rbind(
    ss[which(ss$analysis == paste("trio_concordance", vc, sep=".")
             & ss$measure == "median_child_proportion_inherited"), ],
    ss[which(ss$analysis == paste("trio_concordance", vc, sep=".")
             & ss$measure == "child_proportion_inherited_reciprocal_dynamic_range"), ],
    ss[which(ss$analysis == paste("gt_benchmarking.twin_replicate", vc, sep=".")
             & ss$measure == "median_match_rate"), ],
    ss[which(ss$analysis == paste("gt_benchmarking.twin_replicate", vc, sep=".")
             & ss$measure == "match_rate_reciprocal_dynamic_range"), ],
    ss[which(ss$analysis == paste("heterozygosity", vc, sep=".")
             & ss$measure == "median"), ],
    ss[which(ss$analysis == paste("heterozygosity", vc, sep=".")
             & ss$measure == "reciprocal_dynamic_range"), ],
    ss[which(ss$analysis == paste("variants_per_genome", vc, sep=".")
             & ss$measure == "reciprocal_dynamic_range"), ],
    do.call("rbind", lapply(gb.prefixes, function(gbp){
      rbind(ss[which(ss$analysis %in% paste("gt_benchmarking", gbp, vc, sep=".")
                     & ss$measure == "median_sensitivity"), ],
            ss[which(ss$analysis %in% paste("gt_benchmarking", gbp, vc, sep=".")
                     & ss$measure == "sensitivity_reciprocal_dynamic_range"), ],
            ss[which(ss$analysis %in% paste("gt_benchmarking", gbp, vc, sep=".")
                     & ss$measure == "median_ppv"), ],
            ss[which(ss$analysis %in% paste("gt_benchmarking", gbp, vc, sep=".")
                     & ss$measure == "ppv_reciprocal_dynamic_range"), ])
    }))
  ))

  # Assign universal row names
  rownames(gb.df) <-
    unlist(sapply(apply(gb.df[, c("analysis", "measure")], 1, paste, collapse="."),
                  function(qstr){
                    if(endsWith(qstr, "median_child_proportion_inherited")){
                      return("trio_inh_rate")
                    }else if(endsWith(qstr, "child_proportion_inherited_reciprocal_dynamic_range")){
                      return("trio_inh_rate_rdr")
                    }else if(endsWith(qstr, "median_match_rate")){
                      return("rep_match_rate")
                    }else if(endsWith(qstr, "match_rate_reciprocal_dynamic_range")){
                      return("rep_match_rate_rdr")
                    }else if(qstr == paste("heterozygosity", vc, "median", sep=".")){
                      return("heterozygosity")
                    }else if(qstr == paste("heterozygosity", vc, "reciprocal_dynamic_range", sep=".")){
                      return("heterozygosity_rdr")
                    }else if(qstr == paste("variants_per_genome", vc, "reciprocal_dynamic_range", sep=".")){
                      return("var_per_genome_rdr")
                    }
                    for(gbp in gb.prefixes){
                      if(startsWith(qstr, paste("gt_benchmarking", gbp, sep="."))){
                        if(endsWith(qstr, "sensitivity")){
                          return(paste(gbp, "sens", sep="."))
                        }else if(endsWith(qstr, "sensitivity_reciprocal_dynamic_range")){
                          return(paste(gbp, "sens_rdr", sep="."))
                        }else if(endsWith(qstr, "ppv")){
                          return(paste(gbp, "ppv", sep="."))
                        }else if(endsWith(qstr, "ppv_reciprocal_dynamic_range")){
                          return(paste(gbp, "ppv_rdr", sep="."))
                        }
                      }
                    }
                  }))

  return(gb.df)
}

# Format inter-class benchmarking stats for a single VC
get.inter <- function(ss, vc, vcs){
  # Collect data
  inter.df <- as.data.frame(rbind(
    do.call("rbind", lapply(c("heterozygosity_cor", "variant_count_cor"), function(suf){
      do.call("rbind", lapply(setdiff(vcs, "all"), function(vc2){
        ordered.vcs <- if(vc != vc2){intersect(vcs, c(vc, vc2))}else{rep(vc, 2)}
        inter.prefix <- paste(ordered.vcs[1], "vs", ordered.vcs[2], sep="_")
        vc.n <- max(ss$n[which(ss$analysis == paste("variants_per_genome", vc, sep="."))], na.rm=T)
        if(length(vc.n) == 0){
          vc.n <- NA
        }
        if(vc == vc2){
          data.frame("analysis" = paste(inter.prefix, suf, sep="."),
                     "measure" = "r2", "value" = 1, "n" = vc.n)
        }else{
          ss[which(ss$analysis == paste(inter.prefix, suf, sep=".")), ]
        }
      }))
    }))
  ))

  # Assign universal row names
  rownames(inter.df) <-
    unlist(sapply(apply(inter.df[, c("analysis", "measure")], 1, paste, collapse="."),
                  function(qstr){
                    qstr.prefix <- unlist(strsplit(qstr, split=".", fixed=T))[1]
                    vcs.in.qstr <- unique(unlist(strsplit(qstr.prefix, split="_vs_")))
                    vc2.in.qstr <- setdiff(vcs.in.qstr, vc)
                    out.prefix <- if(length(vc2.in.qstr) > 0){vc2.in.qstr[1]}else{vc}
                    if(endsWith(qstr, "heterozygosity_cor.r2")){
                      return(paste(out.prefix, "heterozyg_cor", sep="."))
                    }else if(endsWith(qstr, "count_cor.r2")){
                      return(paste(out.prefix, "count_cor", sep="."))
                    }
                  }))

  return(inter.df)
}

# Load summary statistics & organize in sub-dataframes for plotting
load.ss <- function(tsv.in, ref.prefix=NULL, gb.prefixes=c()){
  ss <- read.table(tsv.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(ss) <- gsub("#", "", colnames(ss), fixed=T)
  ss$analysis <- gsub("\\.all_variants$", ".all", ss$analysis)
  ss <- invert.dynamic.ranges(ss)

  vcs <- c("all", names(var.class.abbrevs))
  ss.l <- lapply(vcs, function(vc){

    # Counts
    count.df <- get.counts(ss, vc)

    # Site benchmarking
    sb.df <- get.sb(ss, vc, ref.prefix)

    # Genotype benchmarking
    gb.df <- get.gb(ss, vc, gb.prefixes)

    # Inter-class benchmarking
    inter.df <- get.inter(ss, vc, vcs)

    # Structured output
    list("counts" = count.df[which(count.df$n > 0 & !is.na(count.df$n)
                                   & !is.infinite(count.df$value)
                                   & !is.na(count.df$value)), ],
         "site_bench" = sb.df[which((sb.df$n > 0 & !is.na(sb.df$n)
                                     & !is.infinite(sb.df$value)
                                     & !is.na(sb.df$value))
                                    | (vc == "all" & sb.df$analysis == "site_ratios.all")), ],
         "gt_bench" = gb.df[which(gb.df$n > 0 & !is.na(gb.df$n)
                                  & !is.infinite(gb.df$value)
                                  & !is.na(gb.df$value)), ],
         "interclass" = inter.df[which(inter.df$n > 0 & !is.na(inter.df$n)
                                       & !is.infinite(inter.df$value)
                                       & !is.na(inter.df$value)), ])
  })
  names(ss.l) <- vcs
  return(ss.l)
}

# Update the map of QC target values based on user-defined custom overrides
update.targets <- function(target.map, tsv.in=NULL){
  if(!is.null(tsv.in)){
    new.df <- read.table(tsv.in, header=F, sep="\t")
    for(i in 1:nrow(new.df)){
      target.map[as.character(new.df[i, 1])] <- as.numeric(new.df[i, 2])
    }
  }
  return(target.map)
}

# Get target values for a given sumstat dataframe
get.targets <- function(ss.df, target.map, default=1){
  targets <- as.numeric(remap(apply(ss.df[, c("analysis", "measure")],
                                    1, paste, collapse=":"),
                              target.map, default=default))
  return(targets)
}

# Evaluate change in values vs. previous run
compare.prev <- function(ss, prev.ss, vc, g, targets, annotate.targets=FALSE,
                         base.color="black", color.mix=0.1){
  vals <- as.numeric(ss[[vc]][[g]]$value)
  n.vals <- length(vals)
  if(is.null(prev.ss)){
    return(list("prev.vals" = rep(NA, n.vals),
                "rel.delta" = rep(NA, n.vals),
                "color" = rep(NA, n.vals),
                "labels" = rep(NA, n.vals),
                "density" = rep(NA, n.vals)))
  }
  prev.vals <- as.numeric(prev.ss[[vc]][[g]][rownames(ss[[vc]][[g]]), "value"])
  delta <- sapply(1:n.vals, function(i){
    # Report relative change towards target, if optioned
    if(annotate.targets & !is.na(targets[i]) & !is.infinite(targets[i])){
      gap <- abs((vals[i] - targets[i]))
      prev.gap <- abs((prev.vals[i] - targets[i]))
      (prev.gap - gap) / prev.gap

      # Otherwise, report relative change vs. previous value
    }else{
      (vals[i] - prev.vals[i]) / prev.vals[i]
    }
  })

  # Organize various other values
  improved <- as.character(delta > 0.005)
  labels <- paste(remap(improved, c("TRUE" = "+", "FALSE" = ""),
                        default.value=NA),
                  round(100 * delta, 0), "%", sep="")
  labels[which(is.na(delta))] <- NA
  colors <- sapply(delta, function(d){
    if(is.na(d)){
      return(NA)
    }else if(d > 0.005){
      col <- boolean.colors["TRUE"]
    }else if(d < -0.005){
      col <- boolean.colors["FALSE"]
    }else{
      col <- "grey50"
    }
    MixColor(base.color, col, amount1=color.mix)
  })
  density <- as.numeric(remap(improved, c("TRUE" = NA, "FALSE" = 20), default.value=NA))

  # Return as list
  return(list("prev.vals" = prev.vals,
              "rel.delta" = delta,
              "colors" = colors,
              "labels" = labels,
              "density" = density))
}


######################
# Plotting Functions #
######################
# Determine vertical dimensions of plot layout
get.layout <- function(ss, margin.vex=0.075, do.layout=TRUE,
                       parmar=c(0, 0.25, 2, 1)){
  # Get strata to be plotted
  vcs <- names(ss)[which(sapply(ss, function(ss.l){nrow(do.call("rbind", ss.l)) > 0}))]
  groups <- unique(unlist(lapply(ss, names)))
  rows <- lapply(groups, function(g){
    unique(unlist(lapply(ss, function(ss.l){
      rownames(ss.l[[g]])
    })))
  })
  names(rows) <- groups
  group.heights <- lapply(rows, length)
  # Only retain groups with at least one row
  groups <- names(which(unlist(group.heights > 0)))
  rows <- rows[groups]
  group.heights <- group.heights[groups]
  n.groups <- length(groups)
  n.rows <- sum(unlist(group.heights))
  margin.buffer <- margin.vex * (n.rows + n.groups - 1)

  # Set layout
  if(do.layout){
    if("counts" %in% groups & n.groups > 1){
      layout(matrix(1:2, byrow=T, nrow=2),
             heights=c(group.heights[["counts"]] + margin.buffer,
                       sum(unlist(group.heights[setdiff(groups, "counts")])) + n.groups - 2 + margin.buffer))
    }else{
      layout(1)
    }
  }

  return(list("vcs" = vcs, "groups" = groups, "group.heights" = group.heights,
              "rows" = rows, "parmar" = parmar))
}

# Helper sub-function to print label aides within plot.left.labels
left.aligned.labels <- function(add.target=FALSE, add.previous=FALSE,
                                y.start=-0.75, rect.offset=-0.1, bar.w.half=0.3){
  axis(2, at=y.start, labels="RDR = reciprocal\ndynamic range", tick=F, line=-1.25,
       cex.axis=4.5/6, las=2, hadj=0, xpd=T)
  y.inc <- y.start + 1.8
  text(x=0.175, y=y.inc, labels="Current\nvalue", pos=2, cex=4.5/6, xpd=T)

  rect(xleft=0.15, xright=0.25,
       ybottom=y.inc - bar.w.half + rect.offset,
       ytop=y.inc + bar.w.half + rect.offset,
       col="gray75", border=NA, bty="n", xpd=T)
  segments(x0=0.25, x1=0.25,
           y0=y.inc - bar.w.half + rect.offset,
           y1=y.inc + bar.w.half + rect.offset,
           lwd=1, col="gray40", lend="square", xpd=T)
  y.inc <- y.inc + 0.6

  if(add.target){
    rect(xleft=0.15, xright=0.25,
         ybottom=y.inc + 0.5 + (0.4/3),
         ytop=y.inc + 0.5 - (0.4/3),
         border=NA, density=40, bty="n", col=boolean.colors[["FALSE"]])
    text(x=0.23, y=y.inc + 0.5, col=boolean.colors[["FALSE"]],
         cex=4.5/6, labels="Gap", pos=4)
    points(x=0.15, y=y.inc + 0.5, pch=10, srt=45, xpd=T)
    text(x=0.15, y=y.inc + 0.5, cex=4.5/6, labels="Target", pos=2)
    y.inc <- y.inc + 1.35
  }

  if(add.previous){
    text(x=0.275, y=y.inc + 0.5, cex=4.5/6, labels="Change from\nprevious", pos=2, xpd=T)
    rect(xleft=c(0.25, 0.3), xright=c(0.3, 0.35),
         ybottom=y.inc + 0.5 - bar.w.half + rect.offset,
         ytop=y.inc + 0.5 + bar.w.half + rect.offset,
         col=adjustcolor(boolean.colors[c("FALSE", "TRUE")], alpha=0.3),
         border=NA, bty="n")
    rect(xleft=c(0.25, 0.3), xright=c(0.3, 0.35),
         ybottom=y.inc + 0.5 - bar.w.half + rect.offset,
         ytop=y.inc + 0.5 + bar.w.half + rect.offset,
         col=boolean.colors[c("FALSE", "TRUE")],
         density=c(20, NA),
         border=boolean.colors[c("FALSE", "TRUE")])
  }
}

# Plot left margin labels
plot.left.labels <- function(ss, ref.title=NULL, sb_prefixes=NULL, sb_titles=NULL,
                             add.previous=FALSE, add.target=FALSE){
  # Get plot parameters
  params <- get.layout(ss)
  params$parmar[4] <- 0

  # Assign row titles
  title.map <- c("count_all" = "Total variants",
                 "count_rare" = "Rare variants",
                 "count_singleton" = "Singletons",
                 "count_per_genome" = "Variants per genome",
                 "hwe_pass" = "Common Hardy-Weinberg pass rate",
                 "hwe_acc" = "AF-weighted HWE accuracy",
                 "ld" = "Common LD tag rate",
                 "site_ratio" = "Class balance",
                 "site_af_cor" = paste("AF correlation vs.", ref.title),
                 "site_sens" = paste("Common sites found from", ref.title),
                 "site_ppv" = paste("Common sites confirmed by", ref.title),
                 "trio_inh_rate" = "Child inheritance rate",
                 "trio_inh_rate_rdr" = "Inheritance rate RDR",
                 "rep_match_rate" = "Twin concordance rate",
                 "rep_match_rate_rdr" = "Twin concordance RDR",
                 "heterozygosity" = "Heterozygosity per genome",
                 "heterozygosity_rdr" = "Heterozygosity RDR",
                 "var_per_genome_rdr" = "Variants per genome RDR",
                 "snv.heterozyg_cor" = "Heterozygosity cor. vs. SNVs",
                 "indel.heterozyg_cor" = "Heterozygosity cor. vs. indels",
                 "sv.heterozyg_cor" = "Heterozygosity cor. vs. SVs",
                 "snv.count_cor" = "Count correlation vs. SNVs",
                 "indel.count_cor" = "Count correlation vs. indels",
                 "sv.count_cor" = "Count correlation vs. SVs")
  if(!is.null(sb_prefixes)){
    if(is.null(sb_titles)){
      sb_titles <- sb_prefixes
    }else{
      # Be smart about when to use sentence case vs. capslock
      sb_titles <- sapply(sb_titles, function(sbt){
        if(substr(sbt, 2, 2) == tolower(substr(sbt, 2, 2))){
          paste(tolower(substr(sbt, 1, 1)), substr(sbt, 2, nchar(sbt)), sep="")
        }else{
          sbt
        }
      })
    }
    metric.map <- c("sens" = "Per-sample recall",
                    "sens_rdr" = "Sample recall RDR",
                    "ppv" = "Per-sample precision",
                    "ppv_rdr" = "Sample precision RDR")
    for(i in 1:length(sb_prefixes)){
      for(metric in c("sens", "sens_rdr", "ppv", "ppv_rdr")){
        field <- paste(sb_prefixes[i], metric, sep=".")
        title.map[field] <- paste(metric.map[metric], "vs.", sb_titles[i])
      }
    }
  }

  # Assign group titles
  group.titles <- c("counts" = "Counts",
                    "site_bench" = "Site benchmarking",
                    "gt_bench" = "Genotype benchmarking",
                    "interclass" = "Inter-class")

  # Add count labels, if present
  did.left.labels <- FALSE
  if("counts" %in% params$groups){
    prep.plot.area(xlims=c(0, 1), ylims=c(params$group.heights$counts, 0),
                   yaxs="r", parmar=params$parmar)
    sapply(1:params$group.heights$counts, function(ri){
      text(x=1, y=ri-0.5, cex=4.5/6, pos=2, xpd=T,
           labels=remap(params$rows$counts[ri], title.map, default.value=NA))
    })
    staple.bracket(x0=0.55, x1=0.55, y0=0, y1=params$group.heights$counts,
                   staple.len=0.03)
    text(x=0.52, y=params$group.heights$counts/2, font=2, cex=5/6, srt=90,
         labels=shorten.text(group.titles["counts"], cex=5/6, font=2,
                             width=params$group.heights$counts,
                             orientation="vertical"))
    left.aligned.labels(add.previous=add.previous, add.target=add.target)
    did.left.labels <- TRUE
  }

  # Add non-count labels
  noncount.groups <- setdiff(params$groups, "counts")
  if(length(noncount.groups) > 0){
    ymax <- sum(unlist(params$group.heights[noncount.groups])) + length(noncount.groups) - 1
    prep.plot.area(xlims=c(0, 1), ylims=c(ymax+0.25, -0.25), parmar=params$parmar)
    last.y <- 0
    for(g in noncount.groups){
      g.start.y <- last.y
      for(r in params$rows[[g]]){
        text(x=1, y=last.y+0.5, cex=4.5/6, pos=2, xpd=T,
             labels=remap(r, title.map, default.value=NA))
        last.y <- last.y + 1
      }
      staple.bracket(x0=0.05, x1=0.05, y0=g.start.y, y1=last.y, staple.len=0.03)
      text(x=0.02, y=mean(c(g.start.y, last.y)), font=2, cex=5/6, srt=90, xpd=T,
           labels=shorten.text(group.titles[g], cex=5/6, font=2,
                               width=last.y - g.start.y, orientation="vertical"))
      last.y <- last.y + 1
      if(!did.left.labels){
        left.aligned.labels(add.previous=add.previous, add.target=add.target)
      }
    }
  }
}

# Plot summary metrics for a single variant class
plot.ss.bars <- function(ss, vc, annotate.targets=TRUE, prev.ss=NULL,
                         bar.color="gray60", bar.vex=0.8){
  # Get plot parameters
  params <- get.layout(ss)
  bar.w.half <- bar.vex / 2
  if(!is.null(prev.ss)){
    params$parmar[4] <- 2.5
  }

  # Plot counts, if present
  if("counts" %in% params$groups){
    # Gather data
    bar.vals <- log10(as.numeric(ss[[vc]]$counts$value))
    bar.mids <- 1:params$group.heights$counts - 0.5
    targets <- log10(get.targets(ss[[vc]]$counts, target.map))
    prev.dat <- compare.prev(ss, prev.ss, vc, "counts", targets,
                             annotate.targets, base.color=bar.color)
    prev.dat$prev.vals <- log10(prev.dat$prev.vals)

    # Prep plot area
    prep.plot.area(xlims=c(0, ceiling(max(c(bar.vals, targets, prev.dat$prev.vals), na.rm=T))),
                   ylims=c(params$group.heights$counts, 0),
                   yaxs="r", parmar=params$parmar)

    # Top X axis
    count.x.at <- axTicks(3)
    count.x.labs <- sapply(10^count.x.at, clean.numeric.labels)
    clean.axis(3, at=count.x.at, labels=count.x.labs, label.units="count",
               tck=-0.045, label.line=-0.75, title.line=0,
               infinite.positive=TRUE,
               title=if(vc=="all"){"All variants"}else{paste(var.class.names[vc], "s", sep="")})

    # Main bars
    rect(xleft=0, xright=bar.vals,
         ybottom=bar.mids - bar.w.half,
         ytop=bar.mids + bar.w.half,
         col=bar.color, border=bar.color)

    # Previous bars & labels
    if(!is.null(prev.ss)){
      prev.notna.idx <- which(!is.na(prev.dat$prev.vals))
      rect(xleft=bar.vals, xright=prev.dat$prev.vals,
           ybottom=bar.mids - bar.w.half,
           ytop=bar.mids + bar.w.half,
           col=adjustcolor(prev.dat$colors, alpha=0.3),
           border=NA, bty="n")
      rect(xleft=bar.vals, xright=prev.dat$prev.vals,
           ybottom=bar.mids - bar.w.half,
           ytop=bar.mids + bar.w.half,
           density=prev.dat$density[prev.notna.idx],
           col=prev.dat$colors[prev.notna.idx],
           border=prev.dat$colors[prev.notna.idx])
      # segments(x0=prev.dat$prev.vals, x1=prev.dat$prev.vals,
      #          y0=bar.mids - bar.w.half, y1=bar.mids + bar.w.half,
      #          lend="square", col=prev.dat$colors)
      text(x=par("usr")[2], y=bar.mids,
           xpd=T, pos=4, cex=4.5/6, font=3, col=prev.dat$colors,
           labels=prev.dat$labels, offset=0.4)
    }

    # Annotate targets
    if(annotate.targets){
      rect(xleft=bar.vals, xright=targets,
           ybottom=bar.mids - (bar.w.half/3),
           ytop=bar.mids + (bar.w.half/3),
           border=NA, density=40, bty="n", col=boolean.colors[["FALSE"]])
      points(x=targets, y=bar.mids, pch=10, xpd=T)
    }

    # Labels + strong current marks
    # segments(x0=bar.vals, x1=bar.vals,
    #          y0=bar.mids - bar.w.half, y1=bar.mids + bar.w.half,
    #          lend="square", lwd=1, col=MixColor(bar.color, "black"))
    label.widths <- bar.vals
    sapply(1:length(bar.vals), function(x){
      bar.label <- clean.numeric.labels(10^bar.vals[x], min.label.length=2)
      if(is.na(bar.vals[x]) | is.infinite(bar.vals[x])){
        return()
      }
      if(strwidth(bar.label, cex=4.5/6) < label.widths[x]){
        bar.label.x <- par("usr")[1] - (0.02*diff(par("usr")[1:2]))
        label.color <- optimize.label.color(bar.color)
      }else{
        bar.label.x <- bar.vals[x] - (0.02*diff(par("usr")[1:2]))
        label.color <- "black"
      }
      text(x=bar.label.x, y=bar.mids[x]+0.05, cex=4.5/6, col=label.color,
           labels=bar.label, pos=4, xpd=T)
    })
  }

  # Plot all other groups
  noncount.groups <- setdiff(params$groups, "counts")
  if(length(noncount.groups) > 0){
    # Get parameters & prep plot area
    ymax <- sum(unlist(params$group.heights[noncount.groups])) + length(noncount.groups) - 1
    prep.plot.area(xlims=c(0, 1), ylims=c(ymax+0.25, -0.25), parmar=params$parmar)
    last.y <- 0

    # Add top Y axis
    clean.axis(3, label.line=-0.75, tck=-0.02)

    # Plot each group
    for(g in noncount.groups){
      # Gather data
      bar.vals <- as.numeric(ss[[vc]][[g]]$value)
      bar.mids <- 1:params$group.heights[[g]] - 0.5 + last.y
      bar.names <- rownames(ss[[vc]][[g]])
      targets <- get.targets(ss[[vc]][[g]], target.map)
      prev.dat <- compare.prev(ss, prev.ss, vc, g, targets,
                               annotate.targets, base.color=bar.color)

      # Custom annotations for class-specific metrics
      sr.idx <- which(rownames(ss[[vc]][[g]]) == "site_ratio")
      sr.map <- c("snv" = "Ti : Tv", "indel" = "ins : del", "sv" = "(DUP+INS) : DEL")
      if(length(sr.idx) > 0){
        text(x=par("usr")[2], y=bar.mids[sr.idx], pos=2, cex=4.5/6, offset=0,
             col=bar.color, labels=sr.map[vc], xpd=T)
      }

      # Bars
      rect(xleft=0, xright=bar.vals,
           ybottom=bar.mids - bar.w.half,
           ytop=bar.mids + bar.w.half,
           col=bar.color, border=bar.color)

      # Previous bars & labels
      if(!is.null(prev.ss)){
        prev.notna.idx <- which(!is.na(prev.dat$prev.vals))
        rect(xleft=bar.vals[prev.notna.idx],
             xright=prev.dat$prev.vals[prev.notna.idx],
             ybottom=bar.mids[prev.notna.idx] - bar.w.half,
             ytop=bar.mids[prev.notna.idx] + bar.w.half,
             col=adjustcolor(prev.dat$colors[prev.notna.idx], alpha=0.3),
             border=NA, bty="n")
        rect(xleft=bar.vals[prev.notna.idx],
             xright=prev.dat$prev.vals[prev.notna.idx],
             ybottom=bar.mids[prev.notna.idx] - bar.w.half,
             ytop=bar.mids[prev.notna.idx] + bar.w.half,
             density=prev.dat$density[prev.notna.idx],
             col=prev.dat$colors[prev.notna.idx],
             border=prev.dat$colors[prev.notna.idx])
        # segments(x0=prev.dat$prev.vals[prev.notna.idx],
        #          x1=prev.dat$prev.vals[prev.notna.idx],
        #          y0=bar.mids[prev.notna.idx] - bar.w.half,
        #          y1=bar.mids[prev.notna.idx] + bar.w.half,
        #          lend="square", col=prev.dat$colors)
        text(x=par("usr")[2], y=bar.mids,
             xpd=T, pos=4, cex=4.5/6, font=3, col=prev.dat$colors,
             labels=prev.dat$labels, offset=0.4)
      }

      # Targets
      if(annotate.targets){
        rect(xleft=bar.vals, xright=targets,
             ybottom=bar.mids - (bar.w.half/3),
             ytop=bar.mids + (bar.w.half/3),
             border=NA, density=40, bty="n", col=boolean.colors[["FALSE"]])
        points(x=targets, y=bar.mids, pch=10, xpd=T)
      }

      # Labels + strong current marks
      # segments(x0=bar.vals, x1=bar.vals,
      #          y0=bar.mids - bar.w.half, y1=bar.mids + bar.w.half,
      #          lend="square", lwd=1, col=MixColor(bar.color, "black"))
      label.widths <- bar.vals
      sapply(1:length(bar.vals), function(x){
        if(is.na(bar.vals[x]) | is.infinite(bar.vals[x])){
          return()
        }
        if(bar.names[x] %in% pct.metrics){
          bar.label <- paste(round(100 * bar.vals[x], 1), "%", sep="")
        }else{
          bar.label <- round(bar.vals[x], 2)
        }
        if(strwidth(bar.label, cex=4.5/6) < as.numeric(label.widths[x])){
          bar.label.x <- par("usr")[1] - (0.02*diff(par("usr")[1:2]))
          label.color <- optimize.label.color(bar.color)
        }else{
          bar.label.x <- bar.vals[x] - (0.02*diff(par("usr")[1:2]))
          label.color <- "black"
        }
        text(x=bar.label.x, y=bar.mids[x]+0.05, cex=4.5/6, col=label.color,
             labels=bar.label, pos=4, xpd=T)
      })

      # Update y pointer
      last.y <- last.y + params$group.heights[[g]] + 1
    }
  }
}

# Main wrapper to plot all summary metrics
plot.ss <- function(ss, out.prefix, prev.ss=NULL, ref.title=NULL,
                    sb_prefixes=NULL, sb_titles=NULL){
  # Determine subplot dimensions
  # Max aspect ratio: 2:1 w:h
  # Illustrator max width: 14" (so 7" at 50% scale)
  # Powerpoint max dims: 10"w x 5.25"h
  params <- get.layout(ss, do.layout=F)
  pdf.height <- 7
  total.width <- c(10, 12, 14)[length(params$vcs)-1]
  left.width <- 3
  vc.width <- (total.width - left.width) / length(params$vcs)

  # Generate one set of plots with and without targets
  for(do.target in c(TRUE, FALSE)){

    # Generate one set of plots with and without previous (if provided)
    for(do.prev in unique(c(FALSE, !is.null(prev.ss)))){
      # Set consistent filename suffix
      out.suffix <- paste(if(do.target){"w_targets"}else{NULL},
                          if(do.prev){"w_previous"}else{NULL},
                          sep=".")
      out.suffix <- gsub("\\.$", "", gsub("^\\.", "", gsub("[\\.]+", ".", out.suffix)))

      # Plot left axis titles
      pdf(gsub("[\\.]+", ".", paste(out.prefix, "legend", out.suffix, "pdf", sep=".")),
          height=pdf.height, width=left.width)
      plot.left.labels(ss, ref.title, sb_prefixes, sb_titles,
                       add.previous=do.prev, add.target=do.target)
      dev.off()

      sapply(params$vcs, function(vc){
        pdf(gsub("[\\.]+", ".", paste(out.prefix, vc, out.suffix, "pdf", sep=".")),
            height=pdf.height, width=vc.width)
        plot.ss.bars(ss, vc, prev.ss=if(do.prev){prev.ss}else{NULL},
                     bar.color=remap(vc, var.class.colors, default="gray50"),
                     annotate.targets=do.target)
        dev.off()
      })
    }
  }
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot overall QC summary")
parser$add_argument("--stats", metavar=".tsv", type="character",
                    help="File of all QC summary statistics", required=TRUE)
parser$add_argument("--previous-stats", metavar=".tsv", type="character",
                    help=paste("Optional file of QC summary statistics from",
                               "prior run (for comparison)"))
parser$add_argument("--site-ref-prefix", metavar="string", type="character",
                    help="String prefix for site-level benchmarking metrics")
parser$add_argument("--site-ref-title", metavar="string", type="character", default="Ref. cohort",
                    help="Reference dataset title for site-level labels in plots")
parser$add_argument("--sample-benchmarking-prefix", metavar="string",
                    type="character", action="append",
                    help="String prefix for each sample-level benchmarking dataset")
parser$add_argument("--sample-benchmarking-title", metavar="string",
                    type="character", action="append",
                    help="Plot title for each sample-level benchmarking dataset")
parser$add_argument("--custom-targets", metavar=".tsv", type="character",
                    help="Optional two-column .tsv specifying user-defined targets")
parser$add_argument("--custom-constants", metavar=".R", type="character",
                    help="Optional file of custom constants to use for plotting")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV (SINGLE CLASS)
# args <- list("stats" = "~/Downloads/dfci-g2c.v1.gatksv_qc_post_imputation.all_qc_summary_metrics.tsv.gz",
#              "previous_stats" = "~/Downloads/dfci-g2c.v1.initial_gatksv_qc.all_qc_summary_metrics.tsv.gz",
#              "site_ref_prefix" = "gnomad-sv_v4.1",
#              "site_ref_title" = "gnomAD-SV v4.1",
#              "sample_benchmarking_prefix" = c("external_srwgs", "external_lrwgs"),
#              "sample_benchmarking_title" = c("External srWGS", "External lrWGS"),
#              "custom_targets" = "~/scratch/dfci-g2c.v1.gatksv.qc_targets.tsv",
#              "custom_constants" = NULL,
#              "out_prefix" = "~/scratch/dfci-g2c.v1.gatk-sv.initial_qc")

# Load custom constants if optioned
if(!is.null(args$custom_constants)){
  source(args$custom_constants)
}

# Load and organize summary stats
ss <- load.ss(args$stats, args$site_ref_prefix, args$sample_benchmarking_prefix)
if(!is.null(args$previous_stats)){
  prev.ss <- load.ss(args$previous_stats, args$site_ref_prefix, args$sample_benchmarking_prefix)
}else{
  prev.ss <- NULL
}

# Update some constants defined above based on user inputs
if(!is.null(args$sample_benchmarking_prefix)){
  new.pct.metrics <- unlist(sapply(args$sample_benchmarking_prefix,
                                   function(s){paste(s, c("sens", "ppv"), sep=".")}))
  pct.metrics <- unique(c(pct.metrics, new.pct.metrics))
}
target.map <- update.targets(target.map, args$custom_targets)

# Plot summary stats
plot.ss(ss, args$out_prefix, prev.ss=prev.ss, args$site_ref_title,
        args$sample_benchmarking_prefix, args$sample_benchmarking_title)

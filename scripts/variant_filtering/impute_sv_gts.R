#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Impute genotypes for a single SV


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(caret, quietly=TRUE)
require(G2CR, quietly=TRUE)


##################
# Data functions #
##################
# Load & normalize sample covariates
load.covars <- function(tsv.in, keep.samples=NULL){
  if(is.null(tsv.in)){
    return(NULL)
  }
  covars <- read.table(tsv.in, header=T, sep="\t", comment.char="")
  rownames(covars) <- covars[, 1]
  covars[, 1] <- NULL
  if(!is.null(keep.samples)){
    covars <- covars[intersect(rownames(covars), keep.samples), ]
  }
  for(i in 1:ncol(covars)){
    covars[, i] <- scale(as.numeric(covars[, i]))
  }
  impute.missing.values(covars)
}

# Load group information and filter for group-specific imputation eligibility
load.groups <- function(groups.tsv, target.sids, sv.ad, min.ac){
  # Load group info
  if(!is.null(groups.tsv)){
    g.df <- read.table(groups.tsv, header=F, sep="\t")
    g.df <- g.df[which(g.df[, 1] %in% target.sids), ]
    groups <- unique(g.df[, 2])
  }else{
    g.df <- data.frame(target.sids, "ALL")
    groups <- "ALL"
  }
  g.mems <- sapply(groups, function(gid){
    as.character(g.df[which(g.df[, 2] == gid), 1])
  })

  # Count number of ref & alt SV alleles per group
  g.k <- as.data.frame(do.call("rbind", lapply(groups, function(gid){
    c(sum(sv.ad[intersect(names(sv.ad), g.mems[[gid]])] == 0, na.rm=T),
      sum(sv.ad[intersect(names(sv.ad), g.mems[[gid]])] > 0, na.rm=T))
  })))
  colnames(g.k) <- c("ref", "alt")
  rownames(g.k) <- groups
  g.k <- g.k[order(-apply(g.k, 1, min)), ]

  # Iteratively collapse groups until all groups are >= min.ac
  if(any(apply(g.k, 1, min) < min.ac)){
    g.k["remaining", ] <- c(0, 0)
    other.g <- c()
    while(any(apply(g.k, 1, min) < min.ac)){
      next.g <- tail(setdiff(rownames(g.k), "remaining"), 1)
      other.g <- c(other.g, next.g)
      g.k["remaining", ] <- g.k["remaining", ] + g.k[next.g, ]
      g.k <- g.k[setdiff(rownames(g.k), next.g), ]
      g.k <- g.k[order(-apply(g.k, 1, min)), ]
    }
    groups <- rownames(g.k)
    g.mems[["remaining"]] <- unique(unlist(g.mems[other.g]))
    g.mems <- g.mems[groups]
  }

  # Return named list of group memberships
  names(g.mems) <- groups
  return(g.mems)
}

# Create cross-validation folds balanced by non-ref GT count
make.cv.folds <- function(ad.v, k=5, seed=2025){
  set.seed(seed)
  ref.folds <- createFolds(which(ad.v == 0), k)
  alt.folds <- createFolds(which(ad.v > 0), k)
  lapply(1:k, function(i){
    sort(unique(c(ref.folds[[i]], alt.folds[[i]])))
  })
}

# Adjust SV allele dosages for covariates
adjust.sv.ad <- function(sv.ad, covars, train.sids, groups, sv.id, k=5, seed=2025){
  # If covariates aren't supplied, do nothing
  if(is.null(covars)){
    return(sv.ad)
  }

  # Otherwise, train covariate adjustment model on best-powered group
  # We assume technical factors that influence genotype distributions are
  # not causally related to genetic ancestry/group membership
  best.group <- names(groups)[1]
  cov.train.sids <- intersect(train.sids, groups[[1]])
  cat(paste(" - Adjusting raw allele dosages for", sv.id,
            "based on supplied covariates for",
            prettyNum(length(cov.train.sids), big.mark=","),
            best.group, "samples...\n"))
  train.covars <- covars[cov.train.sids, ]
  cov.train.sv.ad <- sv.ad[cov.train.sids]
  cov.train.fold.indexes <- make.cv.folds(cov.train.sv.ad, k=k, seed=seed)
  cov.fit <- train.elastic.net.cv(train.covars, cov.train.sv.ad,
                                  fold.indexes=cov.train.fold.indexes,
                                  seed=seed, tune.length=10)

  # Adjust and return all SV ADs
  cov.fit.e <- mean(coef(cov.fit$finalModel)["(Intercept)", ])
  pred.ad <- predict(cov.fit, covars)

  # Report variance explained
  raw.ad <- sv.ad[intersect(rownames(covars), names(sv.ad))]
  rmse <- round(sqrt(mean((raw.ad - pred.ad)^2)), 3)
  r2 <- cor(raw.ad, pred.ad)^2
  pct.var <- paste(round(100 * r2, 1), "%", sep="")
  cat(paste( " - Covariates explained ",  pct.var,
             " variance in raw SV ADs (RMSE = ", rmse,
             ")\n", sep=""))

  # Return adjusted SV ADs
  adj.ad <- sv.ad[intersect(rownames(covars), names(sv.ad))] - (pred.ad - cov.fit.e)
  return(adj.ad)
}

# Train an SV genotype imputation model
train.imputation <- function(sv.ad, snp.ad, train.sids, k=5, seed=2025){
  # Subset input data to training sids
  train.sv.ad <- sv.ad[intersect(train.sids, names(sv.ad))]
  train.snp.ad <- snp.ad[intersect(train.sids, rownames(snp.ad)), , drop=FALSE]

  # Split samples into cross-validation folds balanced by number of non-ref GTs
  fold.indexes <- make.cv.folds(round(train.sv.ad), k=k, seed=seed)

  # Fit imputation model
  train.elastic.net.cv(train.snp.ad, train.sv.ad, fold.indexes=fold.indexes, seed=seed)
}

# Impute & scale SV allele dosages from a trained allele dosage regression model
impute.sv.ads <- function(sv.fit, snp.ad, train.sv.ad, seed=2025){
  # Apply model to get raw SV ADs
  pred.ad <- predict(sv.fit, snp.ad)

  # Get centroid initialization for AC=0,1,2 (if observed)
  # Note that this is only used for scaling, so doesn't need to be perfect
  # We exclude any het/hom samples with predicted ADs in the bottom 5%,
  # as these are most likely misgenotyped ref samples
  elig.nonref <- names(which(pred.ad > 0.05*diff(range(pred.ad))))
  elig.train.ids <- which(train.sv.ad == 0 | names(train.sv.ad) %in% elig.nonref)
  obs.acs <- intersect(0:2, unique(train.sv.ad[elig.train.ids]))
  train.sv.ac <- sapply(obs.acs, function(k){length(which(train.sv.ad[elig.train.ids] == k))})
  ref.start <- if(0 %in% obs.acs & train.sv.ac[1] > 0){
    median(pred.ad[names(which(train.sv.ad == 0))])
  }else{0}
  het.start <- hom.start <- NULL
  het.start <- if(1 %in% obs.acs){
    if(train.sv.ac[which(obs.acs == 1)] > 0){
      median(pred.ad[intersect(names(which(train.sv.ad == 1)), elig.nonref)])
    }
  }
  if(is.null(het.start)){
    het.start <- sum(range(pred.ad, na.rm=T))/2
  }
  hom.start <- if(2 %in% obs.acs){
    if(train.sv.ac[which(obs.acs == 2)] > 0){
      median(pred.ad[intersect(names(which(train.sv.ad == 2)), elig.nonref)])
    }
  }
  if(is.null(hom.start)){
    hom.start <- max(pred.ad, na.rm=T)
  }
  # If hom.start is not naturally >20% greater than het.start,
  # this likely indicates poor input genotypes that do not clearly distinguish
  # between het and hom clearly. In this case, we pretend no hom GTs exist
  if(hom.start <= 1.2*het.start){
    hom.start <- NULL
    obs.acs <- setdiff(obs.acs, 2)
  }else{
    # Otherwise, ensure hom.start is at least 80% greater than het start
    hom.start <- max(hom.start, het.start + (0.8*(het.start-ref.start)))
  }
  k.start <- c(ref.start, het.start, hom.start)[obs.acs+1]

  # Cluster all samples in untransformed AD space to identify high-confidence samples
  set.seed(seed)
  pred.ac <- obs.acs[kmeans(pred.ad, centers=k.start)$cluster]
  names(pred.ac) <- names(pred.ad)
  k.sids <- sapply(obs.acs, function(ac){
    intersect(names(which(pred.ac == ac)),
              names(which(train.sv.ad == ac)))
  })
  scale.obs.acs <- obs.acs[which(sapply(k.sids, length) > 0)]

  # Scale & return predicted ADs
  if(0 %in% scale.obs.acs){
    pred.ad <- pred.ad - median(pred.ad[k.sids[[1]]])
  }
  if(all(c(1, 2) %in% scale.obs.acs)){
    het.scalar <- 1 / median(pred.ad[k.sids[[2]]])
    hom.scalar <- 2 / median(pred.ad[k.sids[[3]]])
    scalar <- weighted.mean(c(het.scalar, hom.scalar),
                            sqrt(sapply(k.sids[2:3], length)))
  }
  if(2 %in% scale.obs.acs){
    scalar <- 2 / median(pred.ad[k.sids[[3]]])
  }else if(1 %in% scale.obs.acs){
    scalar <- 2 / median(2*pred.ad[k.sids[[2]]])
  }else{
    scalar <- 2
  }
  pred.ad * scalar
}

# Compute P-values for estimated allele dosages versus a parameterized Gaussian
gt.pval <- function(ads, gt.d){
  g.mean <- gt.d[1]
  g.sd <- gt.d[2]
  z <- (ads - g.mean) / g.sd
  (2 * pnorm(abs(z), lower.tail=F))
}

impute.gts <- function(pred.ad, train.sv.ad, min.n.per.ac=10, default.sd=0.2){
  # Define set of samples with concordant imputed and genotyped SV AC
  true.gts <- train.sv.ad[which(train.sv.ad == round(pred.ad[names(train.sv.ad)]))]
  k.sids <- sapply(0:2, function(ac){names(which(true.gts == ac))})

  # Parameterize Gaussians for assigning genotype
  ref.g <- het.g <- hom.g <- NULL
  if(length(k.sids[[1]]) > min.n.per.ac){
    ref.g <- c(mean(pred.ad[k.sids[[1]]]), sd(pred.ad[k.sids[[1]]]))
  }
  if(length(k.sids[[2]]) > min.n.per.ac){
    het.g <- c(mean(pred.ad[k.sids[[2]]]), sd(pred.ad[k.sids[[2]]]))
  }
  if(length(k.sids[[3]]) > min.n.per.ac){
    hom.g <- c(mean(pred.ad[k.sids[[3]]]), sd(pred.ad[k.sids[[3]]]))
  }
  if(is.null(ref.g)){
    if(!is.null(het.g)){
      ref.g <- c(0, het.g[2])
    }else if(!is.null(hom.g)){
      ref.g <- c(0, hom.g[2])
    }else{
      ref.g <- c(0, default.sd)
    }
  }
  if(is.null(het.g)){
    if(!is.null(hom.g)){
      het.g <- c(hom.g[1]/2, hom.g[2])
    }else if(!is.null(ref.g)){
      het.g <- c(1, ref.g[2])
    }else{
      het.g <- c(1, default.sd)
    }
  }
  if(is.null(hom.g)){
    if(!is.null(het.g)){
      hom.g <- c(2*het.g[1], het.g[2])
    }else if(!is.null(ref.g)){
      hom.g <- c(2, ref.g[2])
    }else{
      hom.g <- c(2, default.sd)
    }
  }

  # Assign GT PL to each sample according to GATK formulation
  # See: https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs
  gt.pl <- data.frame("ref" = -10*log10(gt.pval(pred.ad, ref.g)),
                      "het" = -10*log10(gt.pval(pred.ad, het.g)),
                      "hom" = -10*log10(gt.pval(pred.ad, hom.g)))
  pl.v <- as.numeric(unlist(gt.pl))
  pl.max <- 2 * max(pl.v[which(!is.infinite(pl.v))], na.rm=T)
  for(k in 1:3){
    gt.pl[which(is.infinite(gt.pl[, k])), k] <- pl.max
  }

  # Normalize PL per sample
  gt.pl.norm <- as.data.frame(t(apply(gt.pl, 1, function(v){
    v - min(v, na.rm=T)
  })))

  # Compute GT and GQ per sample
  gt.gq <- t(apply(gt.pl.norm, 1, function(pls){
    best.idx <- head(which(pls == min(pls)), 1)
    gt <- c("0/0", "0/1", "1/1")[best.idx]
    gq <- round(min(c(abs(min(pls[-best.idx]) - pls[best.idx]), 99)), 0)
    c(gt, gq)
  }))

  # Summarize imputation results
  res <- as.data.frame(merge(gt.gq, as.data.frame(pred.ad), by="row.names"))
  colnames(res) <- c("sample", "GT", "GQ", "AD")
  res$AD <- round(res$AD, 2)
  return(res)
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Impute GTs for a single SV")
parser$add_argument("--ad", metavar=".tsv", type="character", required=TRUE,
                    help="Allele dosage matrix generated by extract_ad_matrix.py")
parser$add_argument("--sv-id", metavar="string", type="character", required=TRUE,
                    help="Variant ID for the SV to be imputed")
parser$add_argument("--sample-covariates", metavar=".tsv", type="character",
                    help="Optional .tsv of sample covariates for training")
parser$add_argument("--sample-group-labels", metavar=".tsv", type="character",
                    help=paste("Two-column .tsv mapping sample IDs to major",
                               "group labels, like ancestry. If supplied, will",
                               "attempt imputation within each group that meets",
                               "--min-ac requirement"))
parser$add_argument("--min-ac", metavar="int", type="numeric", default=20,
                    help=paste("Minimum number of AC and ref alleles per group",
                               "in --sample-group-labels to permit",
                               "group-specific training."))
parser$add_argument("--min-accuracy", metavar="float", type="numeric", default=0.5,
                    help=paste("Minimum accuacy for carrier status between",
                               "original and imputed genotypes to accept the",
                               "imputation model as trustworthy. Models with",
                               "carrier accuracy below this threshold will only",
                               "write a header to --out-tsv."))
parser$add_argument("--min-r2", metavar="float", type="numeric", default=0.3,
                    help=paste("Minimum coefficient of determination (R2) for ",
                               "adjusted/original and imputed SV allele dosages",
                               "to accept the imputation model as trustworthy.",
                               "Models with R2 below this threshold will only",
                               "write a header to --out-tsv."))
parser$add_argument("--out-tsv", metavar="path", type="character", required=TRUE,
                    help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("ad" = "~/Downloads/dfci-g2c.v1.chr19.final_cleanup_DEL_chr19_4680.ad.tsv.gz",
#              "sv_id" = "dfci-g2c.v1.chr19.final_cleanup_DEL_chr19_4680",
#              "sample_covariates" = "~/Downloads/dfci-g2c.v1.sv_imputation_covariates.tsv.gz",
#              "sample_group_labels" = "~/scratch/dfci-g2c.v1.qc_ancestry.tsv",
#              "min_ac" = 50,
#              "min_accuracy" = 0.7,
#              "min_r2" = 0.3,
#              "out_tsv" = "~/scratch/sv_imp.test.tsv")
# args <- list("ad" = "~/Downloads/dfci-g2c.v1.chr19.final_cleanup_INS_chr19_29.ad.tsv.gz",
#              "sv_id" = "dfci-g2c.v1.chr19.final_cleanup_INS_chr19_29",
#              "sample_covariates" = "~/Downloads/dfci-g2c.v1.sv_imputation_covariates.tsv.gz",
#              "sample_group_labels" = "~/scratch/dfci-g2c.v1.qc_ancestry.tsv",
#              "min_ac" = 50,
#              "min_accuracy" = 0.7,
#              "min_r2" = 0.3,
#              "out_tsv" = "~/scratch/sv_imp.test.tsv")

# Initialize reporting
cat(paste("\n\nNow starting GT imputation for ", args$sv_id, "...\n", sep=""))

# Load allele dosage matrix and split SV from SNPs
cat(paste(" - Loading data for ", args$sv_id, "...\n", sep=""))
ad <- read.table(args$ad, header=T, sep="\t", comment.char="", check.names=F)
rownames(ad) <- ad$sample
ad$sample <- NULL
sv.ad <- ad[, args$sv_id]
names(sv.ad) <- rownames(ad)
snp.ad <- as.data.frame(ad[, setdiff(colnames(ad), args$sv_id)])
rownames(snp.ad) <- rownames(ad)
snp.ad <- impute.missing.values(snp.ad)
target.sids <- rownames(snp.ad)
# For glmnet compatability, the snp.ad matrix must have at least two columns
# If only one tag SNP is identified, we add a dummy second SNP with all AD=0
if(ncol(snp.ad) == 1){
  snp.ad$dummy <- 0
}

# If provided, load sample covariates
sv.samples <- names(sv.ad)[which(!is.na(sv.ad))]
covars <- load.covars(args$sample_covariates, keep.samples=sv.samples)

# Define training samples
train.sids <- if(!is.null(covars)){rownames(covars)}else{sv.samples}

# If provided, load groups and evaluate group-specific strata
groups <- load.groups(args$sample_group_labels, target.sids, sv.ad, args$min_ac)

# Adjust SV ADs if covariates were provided
sv.ad.adj <- adjust.sv.ad(sv.ad, covars, train.sids, groups, args$sv_id)

# Train & apply imputation model for each group
pred.ad <- c()
for(group in names(groups)){
  n.samples <- length(groups[[group]])
  cat(paste(" - Imputing SV allele dosages for",
            prettyNum(n.samples, big.mark=","), group, "samples from",
            prettyNum(length(setdiff(colnames(snp.ad), "dummy")), big.mark=","),
            "tag SNPs...\n"))
  g.train.ids <- intersect(train.sids, groups[[group]])
  group.fit <- train.imputation(sv.ad.adj, snp.ad, g.train.ids)
  g.pred.ad <- impute.sv.ads(group.fit,
                             snp.ad[intersect(rownames(snp.ad), groups[[group]]), ],
                             sv.ad[g.train.ids])
  pred.ad <- c(pred.ad, g.pred.ad)
}

# Probabilistic genotype assignment
cat(paste(" - Predicting SV genotypes for",
          prettyNum(length(target.sids), big.mark=","),
          "samples...\n"))
imp.res <- impute.gts(pred.ad, sv.ad[train.sids])

# Check concordance of imputed GTs vs. original genotypes
c.dat <- merge(data.frame("sv.ad"=sv.ad[which(!is.na(sv.ad))]), imp.res,
               by.x="row.names", by.y="sample",
               all=F, sort=F)[, c("sv.ad", "AD", "GT")]
final.r2 <- cor(c.dat$sv.ad, c.dat$AD, use="complete.obs")^2
c.dat$OGT <- remap(c.dat$sv.ad, c("0" = "0/0", "1" = "0/1", "2" = "1/1"))
gt.acc <- sum(c.dat$OGT == c.dat$GT) / nrow(c.dat)
carrier.acc <- sum(apply(c.dat[, c("GT", "OGT")], 1, function(v){
  length(unique(v == "0/0")) == 1
})) / nrow(c.dat)

# Report confusion matrix for logging
cat(" - Confusion matrix of imputed vs. original GTs:\n")
cm.df <- c.dat[, c("GT", "OGT")]
colnames(cm.df) <- c("Imputed GT", "Original GT")
cat(paste("   ", capture.output(table(cm.df)), "\n", sep=""))

# Write imputed genotype information to --out-tsv
# Only write genotypes if final carrier accuracy is >= --min-accuracy
imp.res$`#sv_id` <- args$sv_id
if(carrier.acc >= args$min_accuracy & final.r2 >= args$min_r2){
  cat(paste(" - Imputation model well-fit; accepting imputed genotypes.\n",
            "   (Carrier accuracy = ", round(carrier.acc, 2),
            "; GT accuracy = ", round(gt.acc, 2),
            "; AD R2 = ", round(final.r2, 2),
            ")\n", sep=""))
}else{
  cat(paste(" - Imputation model below acceptable performance; rejecting imputed genotypes.\n",
            "   (Carrier accuracy = ", round(carrier.acc, 2),
            "; GT accuracy = ", round(gt.acc, 2),
            "; AD R2 = ", round(final.r2, 2),
            ")\n", sep=""))
  imp.res <- imp.res[-(1:nrow(imp.res)), ]
}
write.table(imp.res[, c("#sv_id", setdiff(colnames(imp.res), "#sv_id"))],
            args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)


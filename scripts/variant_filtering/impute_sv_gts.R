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

# Determine the group of sample IDs that has the best balance of ref and nonref SV GTs
find.best.train.group <- function(sv.ad, groups.tsv){
  g.df <- read.table(groups.tsv, header=F, sep="\t")
  groups <- unique(g.df[, 2])
  g.k <- sapply(groups, function(gid){
    g.mems <- as.character(g.df[which(g.df[, 2] == gid), 1])
    min(table(sv.ad[intersect(names(sv.ad), g.mems)] > 0))
  })
  best.g <- head(names(g.k == max(g.k, na.rm=T)), 1)
  as.character(g.df[which(g.df[, 2] == best.g), 1])
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

# Train an SV genotype imputation model
train.imputation <- function(sv.ad, snp.ad, train.sids, covars=NULL,
                             cov.train.sids=train.sids, k=5, seed=2025){
  # Subset input data to training sids
  train.sv.ad <- sv.ad[intersect(train.sids, names(sv.ad))]
  train.snp.ad <- snp.ad[intersect(train.sids, rownames(snp.ad)), ]

  # Split samples into cross-validation folds balanced by number of non-ref GTs
  fold.indexes <- make.cv.folds(train.sv.ad, k=k, seed=seed)

  # Adjust training SV ADs for technical factors, if optioned
  if(!is.null(covars)){
    train.covars <- covars[cov.train.sids, ]
    cov.train.sv.ad <- train.sv.ad[cov.train.sids]
    cov.train.fold.indexes <- make.cv.folds(cov.train.sv.ad, k=k, seed=seed)
    cov.fit <- train.elastic.net.cv(train.covars, cov.train.sv.ad,
                                    fold.indexes=cov.train.fold.indexes,
                                    seed=seed, tune.length=10)
    cov.fit.e <- mean(coef(cov.fit$finalModel)["(Intercept)", ])
    train.sv.ad <- train.sv.ad - (predict(cov.fit, covars) - cov.fit.e)
  }

  # Fit imputation model
  train.elastic.net.cv(train.snp.ad, train.sv.ad, fold.indexes=fold.indexes, seed=seed)
}

# Compute P-values for estimated allele dosages versus a parameterized Gaussian
gt.pval <- function(ads, gt.d){
  g.mean <- gt.d[1]
  g.sd <- gt.d[2]
  z <- (ads - g.mean) / g.sd
  (2 * pnorm(abs(z), lower.tail=F))
}

# Impute SV genotypes from a trained allele dosage regression model
impute.sv.gts <- function(sv.fit, snp.ad, train.sv.ad, min.n.per.gt=10, seed=2025){
  # Apply model to get raw SV ADs
  pred.ad <- predict(sv.fit, snp.ad)

  # Get centroid initialization for AC=0,1,2
  obs.acs <- intersect(0:2, unique(train.sv.ad))
  k.start <- sapply(obs.acs, function(ac){
    median(pred.ad[names(which(train.sv.ad == ac))])
    })

  # Cluster all samples in untransformed AD space to identify high-confidence samples
  set.seed(seed)
  pred.ac <- kmeans(pred.ad, centers=k.start)$cluster - 1
  k.sids <- sapply(obs.acs, function(ac){
    intersect(names(which(pred.ac == ac)),
              names(which(train.sv.ad == ac)))
    })

  # Scale predicted ADs
  if(0 %in% obs.acs){
    pred.ad <- pred.ad - median(pred.ad[k.sids[[1]]])
  }
  if(2 %in% obs.acs){
    pred.ad <- pred.ad * (2 / median(pred.ad[k.sids[[3]]]))
  }else{
    pred.ad <- pred.ad * (2 / median(2*pred.ad[k.sids[[2]]]))
  }

  # Parameterize Gaussians for assigning genotype
  het.g <- c(mean(pred.ad[k.sids[[2]]]), sd(pred.ad[k.sids[[2]]]))
  ref.g <- NULL
  if(0 %in% obs.acs){
    if(length(k.sids[[1]] > min.n.per.gt)){
      ref.g <- c(mean(pred.ad[k.sids[[1]]]), sd(pred.ad[k.sids[[1]]]))
    }
  }
  if(is.null(ref.g)){
    ref.g <- c(0, het.g[2])
  }
  hom.g <- NULL
  if(2 %in% obs.acs){
    if(length(k.sids[[3]]) > min.n.per.gt){
      hom.g <- c(mean(pred.ad[k.sids[[3]]]), sd(pred.ad[k.sids[[3]]]))
    }
  }
  if(is.null(hom.g)){
    hom.g <- c(2, 1) * het.g
  }

  # Assign GT PL to each sample according to GATK formulation
  # See: https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs
  gt.pl <- data.frame("ref" = -10*log(gt.pval(pred.ad, ref.g)),
                      "het" = -10*log(gt.pval(pred.ad, het.g)),
                      "hom" = -10*log(gt.pval(pred.ad, hom.g)))

  # Normalize PL per sample
  gt.pl.norm <- t(apply(gt.pl, 1, function(v){
    v - min(v, na.rm=T)
  }))

  # Compute GT and GQ per sample
  gt.gq <- t(apply(gt.pl.norm, 1, function(pls){
    best.idx <- which(pls == min(pls))
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
                               "group labels, like ancestry. Will only be used ",
                               "for covariate adjustment during model training."))
parser$add_argument("--out-tsv", metavar="path", type="character", required=TRUE,
                    help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("ad" = "~/Downloads/dfci-g2c.v1.chr19.final_cleanup_DEL_chr19_4680.ad.tsv.gz",
#              "sv_id" = "dfci-g2c.v1.chr19.final_cleanup_DEL_chr19_4680",
#              "sample_covariates" = "~/Downloads/dfci-g2c.v1.sv_imputation_covariates.tsv.gz",
#              "sample_group_labels" = "~/scratch/dfci-g2c.v1.qc_ancestry.tsv",
#              "out_tsv" = "~/scratch/sv_imp.test.tsv")

# Load allele dosage matrix and split SV from SNPs
ad <- read.table(args$ad, header=T, sep="\t", comment.char="", check.names=F)
rownames(ad) <- ad$sample
ad$sample <- NULL
sv.ad <- ad[, args$sv_id]
names(sv.ad) <- rownames(ad)
snp.ad <- as.data.frame(ad[, setdiff(colnames(ad), args$sv_id)])
rownames(snp.ad) <- rownames(ad)
snp.ad <- impute.missing.values(snp.ad)

# If provided, load sample covariates
sv.samples <- names(sv.ad)[which(!is.na(sv.ad))]
covars <- load.covars(args$sample_covariates, keep.samples=sv.samples)

# Define training samples
train.sids <- if(!is.null(covars)){rownames(covars)}else{sv.samples}
cov.train.sids <- if(!is.null(covars)){
  intersect(train.sids,
            find.best.train.group(sv.ad[train.sids], args$sample_group_labels))
}else{
  train.sids
}

# Train imputation model
sv.fit <- train.imputation(sv.ad, snp.ad, train.sids, covars, cov.train.sids)

# Apply trained model and predict genotypes for all samples
imp.res <- impute.sv.gts(sv.fit, snp.ad, sv.ad[train.sids])
imp.res$`#sv_id` <- args$sv_id
write.table(imp.res[, c("#sv_id", setdiff(colnames(imp.res), "#sv_id"))],
            args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)


#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Project-wide constants


#' Generate cancer colors
#'
#' Automatically generate colors and palettes for all cancer types in G2C
#'
#' @param cancers Character vector of cancers, in order to be mapped onto rainbow
#' @param n.shades Number of shades (and highlights) appended to each cancer
#' type-specific palette \[default: 2\]
#' @param saturation.range Values passed to [RLCtools::categorical.rainbow()] \[default: c(0.3, 0.95)\]
#' @param value.range Values passed to [RLCtools::categorical.rainbow()] \[default: c(0.35, 1)\]
#' @param period Value passed to [RLCtools::categorical.rainbow()] \[default: `length(cancers) / 3`\]
#' @param plot.colors Should a cancer color diagnostic plot be generated? \[default: FALSE\]
#' @param plot.palettes Should a cancer palette diagnostic plot be generated? \[default: FALSE\]
#'
#' @details A value for "control" should not be provided in `cancers`;
#' this is automatically appended with preset color values
#'
#' @returns two-element list, containing `cancer.colors` and `cancer.palettes`
#'
#' @seealso [RLCtools::hsv.palette()], [RLCtools::categorical.rainbow()]
#'
#' @export create.cancer.colors
#' @export
create.cancer.colors <- function(cancers, n.shades=2, saturation.range=c(0.3, 0.95),
                                 value.range=c(0.35, 1), period=NULL,
                                 plot.colors=FALSE, plot.palettes=FALSE){
  require(RLCtools, quietly=TRUE)

  # Enforce pan-cancer as the first color, always
  cancers <- c("pancan", setdiff(cancers, "pancan"))

  # Get required parameters
  n.cancers <- length(cancers)
  if(is.null(period)){
    period <- n.cancers / 3
  }

  # Generate main cancer color palette
  cancer.colors <- categorical.rainbow(n.cancers, saturation.range=saturation.range,
                                       value.range=value.range, period=period)
  names(cancer.colors) <- cancers
  cancer.colors["oral_cavity"] <- cancer.colors["oral"]
  cancer.colors["control"] <- "#D6D6D6"
  cancer.colors["multiple"] <- cancer.colors["other"] <- cancer.colors["pancan"]
  cancer.colors[c("unknown", "not_specified", "NA")] <- "gray95"
  cancer.colors[c("pancan", "all", "case")] <- "#C43825"

  # Visualize cancer colors to screen, if optioned
  if(plot.colors){
    plot(length(cancer.colors):1, 1:length(cancer.colors), pch=18, cex=4,
         col=rev(cancer.colors), xlim=c(0, length(cancer.colors)+4),
         ylim=c(0, 1.15*length(cancer.colors)), xaxt="n", xlab="", yaxt="n", ylab="")
    text(x=length(cancer.colors):1, y=1:length(cancer.colors), pos=4,
         labels=rev(names(cancer.colors)), srt=45)
  }

  # Generate one palette for each cancer type
  cancer.palettes <- lapply(cancer.colors, function(color){
    pal.w.buffers <- colorRampPalette(c("black", color, "white"))((2*n.shades)+3)
    p <- pal.w.buffers[-c(1, length(pal.w.buffers))]
    p.names <- "main"
    if(n.shades > 0){
      p.names <- c(paste("dark", n.shades:1, sep=""),
                   p.names,
                   paste("light", 1:n.shades, sep=""))
    }
    names(p) <- p.names
    return(p)
  })

  # Visualize palettes to screen, if optioned
  if(plot.palettes){
    RLCtools::prep.plot.area(xlims=c(0, length(cancer.palettes[[1]])+1),
                             ylims=c(0, length(cancer.palettes)),
                             parmar=c(0, 7, 0, 0))
    sapply(1:length(cancer.palettes), function(r){
      points(x=1:length(cancer.palettes[[r]]), y=rep(r, length(cancer.palettes[[r]])),
             pch=19, cex=5, col=cancer.palettes[[r]])
    })
    axis(2, at=(1:length(cancer.palettes))-0.5, las=2, labels=names(cancer.palettes))
  }

  return(list("cancer.colors" = cancer.colors,
              "cancer.palettes" = cancer.palettes))
}


#' Load study constants
#'
#' Load a subset of constants used throughout G2C
#'
#' @param subset Vector of constant groups to load. See `Details` for options. \[default: load all constants\]
#' @param envir Environment passed to [base::assign] \[default: .GlobalEnv\]
#'
#' @details Recognized values for `subset` include:
#' * `colors` : all color palettes used
#' * `scales` : all scales and scale labels
#' * `names` : names of various variables
#' * `other` : miscellaneous other constants
#' * `all` : load all constants
#'
#' @examples
#' # Load list of color palettes
#' get.constants("colors");
#'
#' # Load scales and names colors
#' get.constants(c("scales", "names"))
#'
#' @seealso [base::assign]
#'
#' @export load.constants
#' @export
load.constants <- function(subset="all", envir=.GlobalEnv){
  # Define colors
  cancer.color.order <- c("pancan", "sarcoma", "oral", "melanoma", "esophagus",
                          "thyroid","lung", "liver", "kidney", "bladder", "cns",
                          "colorectal", "prostate", "stomach", "pancreas",
                          "uterus", "ovary", "breast")
  cancer.color.set <- create.cancer.colors(cancer.color.order, period=3.9)
  male.colors <- c("dark2" = "#2A5869",
                   "dark1" = "#3F839D",
                   "main" = "#54AFD1",
                   "light1" = "#76BFDA",
                   "light2" = "#98CFE3")
  female.colors <- c("dark2" = "#73475E",
                     "dark1" = "#AD6B8C",
                     "main" = "#E68EBB",
                     "light1" = "#EBA5C9",
                     "light2" = "#F0BBD6")
  AFR.colors <- c("dark3" = "#796624",
                  "dark2" = "#A18730",
                  "dark1" = "#C9A93C",
                  "main" = "#F1CB48",
                  "light1" = "#F3D160",
                  "light2" = "#F5D878",
                  "light3" = "#F7E091")
  AMR.colors <- c("dark2" = "#652223",
                  "dark1" = "#973234",
                  "main" = "#C94345",
                  "light1" = "#D4696A",
                  "light2" = "#DF8E8F")
  EAS.colors <- c("dark2" = "#5B6519",
                  "dark1" = "#889826",
                  "main" = "#B5CA32",
                  "light1" = "#C4D55B",
                  "light2" = "#D3DF84")
  EUR.colors <- c("dark2" = "#4E6774",
                  "dark1" = "#749AAE",
                  "main" = "#9BCDE8",
                  "light1" = "#AFD7ED",
                  "light2" = "#C3E1F1")
  SAS.colors <- c("dark2" = "#4D2D4E",
                  "dark1" = "#734474",
                  "main" = "#995A9B",
                  "light1" = "#AD7BAF",
                  "light2" = "#C29CC3")
  csq.colors <- c("synonymous" = "#AAAAAA",
                  "missense" = "#FF6103",
                  "lof" = "#9D1309")
  DFCI.colors <- c("darkblue" = "#003354",
                   "midblue" = "#00629B",
                   "lightblue" = "#41B6E6",
                   "paleblue" = "#D4EEFF",
                   "yellow" = "#FFA300")
  var.class.colors <- c("snv" = "#D9C9AE",
                        "indel" = "#A69C88",
                        "sv" = "#736D61")
  snv.colors <- c("ti" = "#BCE3A8",
                  "tv" = "#8ABD71")
  indel.colors <- c("ins" = "#65A1C7",
                    "del" = "#CC6C62")
  sv.colors <- c("DEL" = "#AD574C",
                 "DUP" = "#367FAD",
                 "CNV" = "#7D66AD",
                 "INS" = "#C272A9",
                 "INV" = "#AD8757",
                 "CPX" = "#71B086",
                 "BND" = "#5C7A47",
                 "CTX" = "#5C7A47",
                 "OTH" = "#5C7A47")
  colors <- list(
    "cancer.colors" = cancer.color.set$cancer.colors,
    "cancer.palettes" = cancer.color.set$cancer.palettes,
    "pop.colors" = c("AFR" = AFR.colors[["main"]],
                     "AMR" = AMR.colors[["main"]],
                     "EAS" = EAS.colors[["main"]],
                     "EUR" = EUR.colors[["main"]],
                     "SAS" = SAS.colors[["main"]],
                     "OTH" = "gray"),
    "sex.colors" = c("male" = male.colors[["main"]],
                     "female" = female.colors[["main"]],
                     "other" = "#666245"),
    "male.colors" = male.colors,
    "female.colors" = female.colors,
    "pop.colors" = c("AFR" = AFR.colors[["main"]],
                     "AMR" = AMR.colors[["main"]],
                     "EAS" = EAS.colors[["main"]],
                     "EUR" = EUR.colors[["main"]],
                     "SAS" = SAS.colors[["main"]],
                     "OTH" = "gray"),
    "pop.palettes" = list("AFR" = AFR.colors,
                          "AMR" = AMR.colors,
                          "EAS" = EAS.colors,
                          "EUR" = EUR.colors,
                          "SAS" = SAS.colors,
                          "OTH" = "gray"),
    "AFR.colors" = AFR.colors,
    "AMR.colors" = AMR.colors,
    "EAS.colors" = EAS.colors,
    "EUR.colors" = EUR.colors,
    "SAS.colors" = SAS.colors,
    "subpop.colors" = c("ACB" = AFR.colors[["light1"]],
                        "ASW" = AFR.colors[["main"]],
                        "ESN" = AFR.colors[["dark3"]],
                        "GWD" = AFR.colors[["light2"]],
                        "LWK" = AFR.colors[["dark2"]],
                        "MSL" = AFR.colors[["light3"]],
                        "YRI" = AFR.colors[["dark1"]],
                        "CLM" = AMR.colors[["light2"]],
                        "MXL" = AMR.colors[["light1"]],
                        "PEL" = AMR.colors[["dark1"]],
                        "PUR" = AMR.colors[["main"]],
                        "CDX" = EAS.colors[["light2"]],
                        "CHB" = EAS.colors[["dark2"]],
                        "CHS" = EAS.colors[["main"]],
                        "JPT" = EAS.colors[["dark1"]],
                        "KHV" = EAS.colors[["light1"]],
                        "CEU" = EUR.colors[["main"]],
                        "FIN" = EUR.colors[["light2"]],
                        "GBR" = EUR.colors[["dark2"]],
                        "IBS" = EUR.colors[["dark1"]],
                        "TSI" = EUR.colors[["light1"]],
                        "BEB" = SAS.colors[["dark1"]],
                        "GIH" = SAS.colors[["light1"]],
                        "ITU" = SAS.colors[["dark2"]],
                        "PJL" = SAS.colors[["main"]],
                        "STU" = SAS.colors[["light2"]]),
    "var.class.colors" = var.class.colors,
    "var.subclass.colors" = c(snv.colors, indel.colors, sv.colors),
    "var.ref.color" = "#99ADBA",
    "snv.colors" = snv.colors,
    "indel.colors" = indel.colors,
    "sv.colors" = sv.colors,
    "csq.colors" = csq.colors,
    "DFCI.colors" = DFCI.colors,
    "stage.colors" = c("0" = "white",
                       "I" = "#F8FAA7",
                       "II" = "#FFCC66",
                       "III" = "#FE8002",
                       "IV" = "#ED3823",
                       "unknown" = "gray85"),
    "relative.colors" = c("duplicates" = "#9D1309",
                          "parent-child" = "#FF6103",
                          "siblings" = "#FFB14D",
                          "unrelated" = "#AAAAAA"),
    "annotation.color" = "gray75",
    "boolean.colors" = c("TRUE" = "#59D05B",
                         "FALSE" = "#f76f57"))

  # Define scales
  logscale.major <- 10^(-10:10)
  contig.lengths <- c("chr1" = 248956422,
                      "chr2" = 242193529,
                      "chr3" = 198295559,
                      "chr4" = 190214555,
                      "chr5" = 181538259,
                      "chr6" = 170805979,
                      "chr7" = 159345973,
                      "chr8" = 145138636,
                      "chr9" = 138394717,
                      "chr10" = 133797422,
                      "chr11" = 135086622,
                      "chr12" = 133275309,
                      "chr13" = 114364328,
                      "chr14" = 107043718,
                      "chr15" = 101991189,
                      "chr16" = 90338345,
                      "chr17" = 83257441,
                      "chr18" = 80373285,
                      "chr19" = 58617616,
                      "chr20" = 64444167,
                      "chr21" = 46709983,
                      "chr22" = 50818468,
                      "chrX" = 156040895,
                      "chrY" = 57227415)
  scales <- list(
    "logscale.major" = logscale.major,
    "logscale.major.bp" = 10^(0:9),
    "logscale.major.bp.labels" = c(sapply(c("bp", "kb", "Mb"),
                                          function(suf){paste(c(1, 10, 100), suf, sep="")}),
                                   "1 Gb"),
    "logscale.demi" = as.numeric(sapply(logscale.major, function(e){c(1, 5)*e})),
    "logscale.demi.bp" = as.numeric(sapply(10^(0:9), function(e){c(1, 5)*e})),
    "logscale.demi.bp.labels" = c(paste(c(1, 5, 10, 50, 100, 500), "bp", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "kb", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "Mb", sep=""),
                                  paste(c(1, 5), "Gb", sep="")),
    "logscale.minor" = as.numeric(sapply(logscale.major, function(e){(1:9)*e})),
    "yearscale.major" = 0:100 * 365,
    "yearscale.demi" = seq(0, 100, 0.5) * 365,
    "yearscale.minor" = seq(0, 100, 1/12) * 365,
    "contig.lengths" = contig.lengths
  )

  # Define names
  all.names <- list(
    "pop.abbreviations" = c("AFR" = "Afr.",
                            "AMR" = "Amer.",
                            "EAS" = "E. Asn.",
                            "EUR" = "Eur.",
                            "SAS" = "S. Asn.",
                            "OTH" = "Other"),
    "pop.names.short" = c("AFR" = "African",
                          "AMR" = "American",
                          "EAS" = "E. Asian",
                          "EUR" = "European",
                          "SAS" = "S. Asian",
                          "OTH" = "Other"),
    "pop.names.long" = c("AFR" = "African/African-American",
                         "AMR" = "Latino/admixed American",
                         "EAS" = "East Asian",
                         "EUR" = "European",
                         "SAS" = "South Asian",
                         "OTH" = "Other/unknown"),
    "sex.names" = c("male" = "Male",
                    "female" = "Female",
                    "other" = "Other"),
    "genetic.sex.names" = c("male" = "XY",
                            "female" = "XX",
                            "other" = "Other"),
    "stage.names" = c("I" = "I",
                      "II" = "II",
                      "III" = "III",
                      "IV" = "IV",
                      "unknown" = "Unknown"),
    "stage.names.long" = c("I" = "Stage I",
                           "II" = "Stage II",
                           "III" = "Stage III",
                           "IV" = "Stage IV",
                           "unknown" = "Unknown"),
    "relative.names" = c("duplicates" = "Identical",
                         "parent-child" = "Parent-child",
                         "siblings" = "Siblings",
                         "unrelated" = "Unrelated"),
    "var.class.names" = c("snv" = "Single nucleotide variant",
                          "indel" = "Small insertion/deletion",
                          "sv" = "Structural variant"),
    "var.class.abbrevs" = c("snv" = "SNV",
                            "indel" = "Indel",
                            "sv" = "SV"),
    "var.subclass.names" = c("ti" = "Transition",
                             "tv" = "Transversion",
                             "del" = "Deletion",
                             "ins" = "Insertion",
                             "DEL" = "Deletion",
                             "DUP" = "Duplication",
                             "INS" = "Insertion",
                             "CNV" = "Multiallelic CNV",
                             "INV" = "Inversion",
                             "CPX" = "Complex SV",
                             "CTX" = "Translocation",
                             "OTH" = "Other SV"),
    "var.subclass.names.long" = c("ti" = "Transition SNV",
                                  "tv" = "Transversion SNV",
                                  "del" = "Small deletion",
                                  "ins" = "Small insertion",
                                  "DEL" = "Large deletion",
                                  "DUP" = "Duplication",
                                  "INS" = "Large insertion",
                                  "CNV" = "Multiallelic copy-number variant",
                                  "INV" = "Inversion",
                                  "CPX" = "Complex structural variant",
                                  "CTX" = "Chromosomal translocation",
                                  "OTH" = "Other structural variant"),
    "var.subclass.names.short" = c("ti" = "Transition",
                                   "tv" = "Transversion",
                                   "del" = "Deletion",
                                   "ins" = "Insertion",
                                   "DEL" = "Deletion",
                                   "DUP" = "Duplication",
                                   "INS" = "Insertion",
                                   "CNV" = "mCNV",
                                   "INV" = "Inversion",
                                   "CPX" = "Complex",
                                   "CTX" = "Transloc.",
                                   "OTH" = "Other"),
    "var.subclass.abbrevs" = c("ti" = "Ti.",
                               "tv" = "Tv.",
                               "del" = "Del.",
                               "ins" = "Ins.",
                               "DEL" = "Del.",
                               "DUP" = "Dup.",
                               "INS" = "Ins.",
                               "CNV" = "mCNV",
                               "INV" = "Inv.",
                               "CPX" = "Cpx.",
                               "CTX" = "Tloc.",
                               "OTH" = "Other"),
    "sv.abbreviations" = c("DEL" = "Del.",
                           "DUP" = "Dup.",
                           "CNV" = "mCNV",
                           "INS" = "Ins.",
                           "INV" = "Inv.",
                           "CPX" = "Cpx.",
                           "CTX" = "Tloc.",
                           "OTH" = "Other"),
    "sv.names" = c("DEL" = "Deletion",
                   "DUP" = "Duplication",
                   "CNV" = "mCNV",
                   "INS" = "Insertion",
                   "INV" = "Inversion",
                   "CPX" = "Complex SV",
                   "CTX" = "Translocation",
                   "OTH" = "Other SV"),
    "csq.names.short" = c("synonymous" = "Syn.",
                          "missense" = "Mis.",
                          "lof" = "LoF"),
    "cancer.names" = c("pancan" = "Pan-cancer",
                       "all" = "All cancer",
                       "case" = "Cancer",
                       "prostate" = "Prostate",
                       "breast" = "Breast",
                       "lung" = "Lung",
                       "colorectal" = "Colorectal",
                       "melanoma" = "Melanoma",
                       "uterus" = "Uterine",
                       "kidney" = "Renal",
                       "bladder" = "Bladder",
                       "oral_cavity" = "Oral",
                       "oral" = "Oral",
                       "ovary" = "Ovarian",
                       "cns" = "CNS",
                       "pancreas" = "Pancreatic",
                       "esophagus" = "Esophageal",
                       "liver" = "Liver",
                       "stomach" = "Gastric",
                       "thyroid" = "Thyroid",
                       "sarcoma" = "Sarcoma",
                       "other" = "Other",
                       "multiple" = "Multiple cancers",
                       "control" = "Control",
                       "unknown" = "Unknown",
                       "not_specified" = "Unknown",
                       "NA" = "Unknown"),
    "cohort.names.vlong" = c("apollo" = "Applied Proteogenomics Organizational Learning and Outcomes",
                             "aou" = "NIH All of Us",
                             "biome" = "BioMe Biobank at Icahn School of Medicine",
                             "ceph" = "Utah Centre d′Etudes du Polymorphisme Humain",
                             "copdgene" = "COPD Genetic Epidemiology",
                             "cptac" = "Clinical Proteomic Tumor Analysis Consortium",
                             "eagle" = "Environment And Genetics in Lung Cancer Etiology",
                             "gmkf" = "Gabriella Miller Kids First",
                             "gtex" = "Genotype-Tissue Expression Project",
                             "hcmi" = "Human Cancer Models Initiative",
                             "hgsvc" = "1000 Genomes Project",
                             "hmf" = "Hartwig Medical Foundation",
                             "icgc" = "International Cancer Genome Consortium",
                             "lcins" = "Zhang et al., Nat. Genet. (2021)",
                             "mesa" = "Multi-Ethnic Study of Atherosclerosis",
                             "proactive-core" = "Dana-Farber Proactive",
                             "proactive-other" = "Dana-Farber other",
                             "proactive" = "Dana-Farber Proactive",
                             "stjude" = "St. Jude Children's Research Hospital",
                             "ufc" = "Dana-Farber Unexplained Familial Cancers",
                             "wcdt" = "NCI West-Coast Dream Team",
                             "other" = "Other cohorts"),
    "cohort.names.long" = c("apollo" = "NCI APOLLO",
                            "aou" = "NIH All of Us",
                            "biome" = "BioMe Biobank",
                            "ceph" = "CEPH families",
                            "copdgene" = "COPDGene",
                            "cptac" = "NCI CPTAC",
                            "eagle" = "NCI EAGLE",
                            "gmkf" = "GMKF",
                            "gtex" = "GTEx",
                            "hcmi" = "HCMI",
                            "hgsvc" = "1000 Genomes Project",
                            "hmf" = "Hartwig Medical Foundation",
                            "icgc" = "ICGC",
                            "lcins" = "Zhang 2021",
                            "mesa" = "MESA",
                            "proactive-core" = "DFCI Proactive",
                            "proactive-other" = "DFCI other",
                            "proactive" = "DFCI Proactive",
                            "stjude" = "St. Jude",
                            "ufc" = "DFCI families",
                            "wcdt" = "NCI WCDT",
                            "other" = "Other"),
    "cohort.names.short" = c("apollo" = "APOLLO",
                             "aou" = "All of Us",
                             "biome" = "BioMe",
                             "ceph" = "CEPH",
                             "copdgene" = "COPDGene",
                             "cptac" = "CPTAC",
                             "eagle" = "EAGLE",
                             "gmkf" = "GMKF",
                             "gtex" = "GTEx",
                             "hcmi" = "HCMI",
                             "hgsvc" = "1000G",
                             "hmf" = "Hartwig",
                             "icgc" = "ICGC",
                             "lcins" = "LCINS",
                             "mesa" = "MESA",
                             "proactive-core" = "Proactive",
                             "proactive-other" = "DFCI other",
                             "proactive" = "Proactive",
                             "stjude" = "St. Jude",
                             "ufc" = "DFCI families",
                             "wcdt" = "WCDT",
                             "other" = "Other"),
    "cohort.type.names" = c("aou" = "All of Us",
                            "cancer" = "Oncology research",
                            "popgen" = "Population genomics"),
    "cohort.type.names.short" = c("aou" = "All of Us",
                                  "cancer" = "Oncology",
                                  "popgen" = "Pop. gen.")
  )

  # Define other constants
  other <- list(
    "cohort.type.map" = c("apollo" = "cancer",
                          "aou" = "aou",
                          "biome" = "popgen",
                          "ceph" = "popgen",
                          "copdgene" = "popgen",
                          "cptac" = "cancer",
                          "eagle" = "cancer",
                          "gmkf" = "cancer",
                          "gtex" = "popgen",
                          "hcmi" = "cancer",
                          "hgsvc" = "popgen",
                          "hmf" = "cancer",
                          "icgc" = "cancer",
                          "lcins" = "cancer",
                          "mesa" = "popgen",
                          "proactive-core" = "cancer",
                          "proactive-other" = "cancer",
                          "proactive" = "cancer",
                          "stjude" = "cancer",
                          "ufc" = "cancer",
                          "wcdt" = "cancer",
                          "other" = "other"),
    "cohort.realigned" = c("apollo" = FALSE,
                           "aou" = FALSE,
                           "biome" = FALSE,
                           "ceph" = TRUE,
                           "copdgene" = FALSE,
                           "cptac" = FALSE,
                           "eagle" = FALSE,
                           "gmkf" = FALSE,
                           "gtex" = FALSE,
                           "hcmi" = FALSE,
                           "hgsvc" = FALSE,
                           "hmf" = TRUE,
                           "icgc" = FALSE,
                           "lcins" = TRUE,
                           "mesa" = FALSE,
                           "proactive-core" = FALSE,
                           "proactive-other" = FALSE,
                           "proactive" = FALSE,
                           "stjude" = FALSE,
                           "ufc" = TRUE,
                           "wcdt" = FALSE)
  )

  # Assign constants to global environment
  if(length(intersect(subset, c("colors", "all"))) > 0){
    for(variable in names(colors)){
      assign(variable, colors[[variable]], envir=envir)
    }
  }
  if(length(intersect(subset, c("scales", "all"))) > 0){
    for(variable in names(scales)){
      assign(variable, scales[[variable]], envir=envir)
    }
  }
  if(length(intersect(subset, c("names", "all"))) > 0){
    for(variable in names(all.names)){
      assign(variable, all.names[[variable]], envir=envir)
    }
  }
  if(length(intersect(subset, c("other", "all"))) > 0){
    for(variable in names(other)){
      assign(variable, other[[variable]], envir=envir)
    }
  }
}

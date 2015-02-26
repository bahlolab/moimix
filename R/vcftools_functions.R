# vcftools_functions.R
# Author: Stuart Lee
# Date: 25-02-2015
# Description: functions for manipulating vcf files
# using the SeqArray package

#' Missing rate per sample
#'
#' @description Calculate the proportion of "No calls" in a sample
#' @details Uses \code{seqApply} to compute the per sample missingness rate
#' from a SeqVarGDSClass object
#' @importFrom SeqArray seqApply
#' @importClassesFrom SeqArray SeqVarGDSClass
#' @param gdsfile a SeqVarGDSClass object
#' @return a numeric vector containing proporition of missing calls in each
#' sample
#' @export
sampleMissingness <- function(gdsfile) {
  stopifnot(class(gdsfile) == "SeqVarGDSClass")
  geno.summary <- seqSummary(gdsfile, "genotype")
  dm <- geno$seldim
  ploidy <- geno$dim[1]
  n <- integer(dm[1])
  seqApply(gdsfile, "genotype",
           function(x) { n <<- n + colSums(is.na(x)) },
           margin = "by.variant", as.is = "none")
  n / (ploidy * dm[2])
}

#' Allele frequency estimate
#'
#' @description Estimate the allele frequency for each variant site
#' @details Uses \code{seqApply} to estimate variant frequency over the
#' SeqVarGDS Class object
#' @importFrom SeqArray seqApply
#' @importClassesFrom SeqArray SeqVarGDSClass
#' @param gdsfile a SeqVarGDSClass object
#' @param reference logical indicating whether to compute the A or B allele
#' frequency (default TRUE)
#' @return a numeric vector containing allele frequency estimates
#' @export
alleleFreqEst <- function(gdsfile, reference = TRUE) {
  stopifnot(class(gdsfile) == "SeqVarGDSClass")

  if(reference) {
    gt <- 0
  } else {
    gt <- 1
  }
  seqApply(gdsfile, "genotype",
           function(x) { mean(x == gt, na.rm = TRUE) },
           margin = "by.variant", as.is = "double")

}

#' P.falciparum chromosome labeller
#' @description Key-value pairs for mapping Pf3D7 chromsome names to numeric
#' chromsome. Used as labeller for facet grid
#' @param variable Pf3D7 chromosome name
#' @param value numeric chromsome name
#' @export
chrom_labeller <- function(variable, value) {

  chrom.names <- list("Pf3D7_01_v3" = "01",
                      "Pf3D7_02_v3" = "02",
                      "Pf3D7_03_v3" = "03",
                      "Pf3D7_04_v3" = "04",
                      "Pf3D7_05_v3" = "05",
                      "Pf3D7_06_v3" = "06",
                      "Pf3D7_07_v3" = "07",
                      "Pf3D7_08_v3" = "08",
                      "Pf3D7_09_v3" = "09",
                      "Pf3D7_10_v3" = "10",
                      "Pf3D7_11_v3" = "11",
                      "Pf3D7_12_v3" = "12",
                      "Pf3D7_13_v3" = "13",
                      "Pf3D7_14_v3" = "14")

  return(chrom.names[value])
}


#' Within sample B-allele frequency estimate
#'
#'@description Estimate within sample B-allele frequency from variant support
#'reads and read depth.
#'@import GenomicRanges
#'@importClassesFrom SeqArray SeqVarGDSClass
#'@importFrom SeqArray seqSetFilter seqGetData
#'@param gdsfile a SeqVarGDS class object
#'@param sample.id character string for sample identifier
#'@param minDP numeric minimum depth to consider
#'@param minGQ numeric minimum geontype quality to consider
#'@param metadata.granges a GRanges class object pertaining to variant sites
#'@return a GRanges class object with metadata slot 'bfreq'
#'@export
bVarFreqEst <- function(gdsfile, sample.id, minDP, minGQ, metadata.granges) {
  # estimate the variant allele frequency on a per sample basis
  # based on coverage, and apply filters
  # reduce gds file to sample of interest
  seqSetFilter(gdsfile, sample.id)
  dp <- seqGetData(gdsfile, "annotation/format/DP")$data
  gq <- seqGetData(gdsfile, "annotation/format/GQ")$data
  ad <- seqGetData(gdsfile, "annotation/format/AD")$data
  filter <- intersect(which(dp > minDP), which(gq > minGQ))
  # create a new Granges Object
  filtered.granges <- metadata.granges[filter]
  elementMetadata(filtered.granges)$bfreq <- ad[filter] / dp[filter]
  # reset filter
  filtered.granges
}

#' Plot B-allele frequency genome wide
#'
#' @description Plot allele-frequency spectrum accross genome and facet by
#' chromosome
#' @importFrom ggbio ggplot
#' @import ggplot2
#' @param sample.granges a GRanges object
#' @return a ggplot2 object
#' @export
plotVarFreq <- function(sample.granges) {
  ggbio::ggplot(sample.granges) +
    geom_point(aes(x = start, y = bfreq), size = 0.5, alpha = 0.05) +
    facet_grid(seqnames ~ ., labeller = chrom_labeller) +
    theme_bw() +
    labs(y = "B-allele frequency",
         x = "Position (bp)") +
    theme(strip.text = element_text(size = 4),
          axis.text.y = element_text(size = 4),
          axis.text.x = element_text(size = 4),
          axis.title.x = element_text(size = 4),
          axis.title.y = element_text(size = 4))
}

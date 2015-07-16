# Title: utility_functions.R
# Author: Stuart Lee, Katherine Smith
# Description: Utility functions for the haldane package
#--- Datasets
#' Minor Allele Frequenies of 86,000 round Plasmodium falciparum SNPs.
#'
#' A dataset containing information about exonic SNPs Plasmodium falciparum
#' from Manske et al., 2012 \url{http://www.ncbi.nlm.nih.gov/pubmed/22722859}.
#'
#' @name global.allele.freq
#' @docType data
#' @format A data frame with 86158 rows and 11 variables:
#' \describe{
#'   \item{ID}{Row number}
#'   \item{Chromosome}{Chromosome}
#'   \item{Position}{Physical position of SNP}
#'   \item{snpName}{SNP identifier}
#'   \item{Geme}{Gene name}
#'   \item{MAF_AFR}{Minor allele frequency in Africa}
#'   \item{MAF_SEA}{Minor allele frequency in South East Asia}
#'   \item{MAF_PNG}{Minor allele frequency in Papua New Guinea}
#'   \item{GeneID}{Gene Identifier}
#'   \item{GeneAliases}{Alternate Gene Identifier}
#'   \item{GeneDescription}{Proposed function of gene}
#'
#'
#' }
#' @source \url{http://www.malariagen.net/data}
NULL

#' Minor allele frequencies of 110,000 PNG Plasmodium falciparum SNPs
#'
#' A dataset containing SNP frequencies derived from 104 isolates
#' from Papua New Guinea
#' @docType data
#' @name png.allele.frq
#' @format A data.frame with 107879 rows and 5 columns
#' \describe{
#'  \item{chr}{Pf3D7 chromsome name}
#'  \item{position}{Physical position}
#'  \item{a.frq}{Reference allele frequency}
#'  \item{b.frq}{Alternate allele frequency}
#'  \item{maf}{Minor allele frequency}
#' }
NULL

# #--- Stastical functions
# #' Harmonic Mean
# #'
# #' This function returns the harmonic mean of an array.
# #' Code from http://economistatlarge.com/r-guide/arithmetic-harmonic-geometric-means-r
# #' @param array numeric vector
# #' @keywords harmonic mean
# #' @export
# #' @examples
# #' harmonicMean(rep(100,10))
# harmonicMean <- function(array) {
# 
# 
#   if (!is.numeric(array)) {
#     stop("Passed argument must be an array. Consider using sapply for data frames.")
#   }
# 
#   if (any(array<=0)) {
#     stop("All values must be greater than zero.")
#   }
#   length(array) / sum(1 / array)
# }
# 
# dirichletAlpha <- function(K,m,a0) {
#   #' Concentration parameters for Dirichlet distribution
#   #'
#   #' @description  Generate vector alpha giving the parameters of a Dirichlet distribution
#   #' @details A starting value a0 is multiplied by a factor m to obtain K categories
#   #' @param K number of categories
#   #' @param m weighting factor
#   #' @param a0 initial value
#   #' @export
#   #'
#   a0*m^(0:(K-1))
# }
# 
# 
# #' Dirichlet distribution random number generator
# #'
# #' @description Draw random numbers from a Dirichlet distribution with parameters given by the vector a.
# #' Code modified from http://stats.stackexchange.com/a/30972
# #'
# #' @param a vector of concentration parameters
# #' @return length(a) dirichlet random numbers
# #' @keywords Dirichlet distribution
# #' @export
# #' @examples
# #' rdirichlet(c(1,2,4))
# rdirichlet <- function(a) {
# 
# 
#   if (!is.numeric(a)) {
#     stop("Argument must be a numeric vector.")
#   }
# 
#   if (any(a<=0)) {
#     stop("All values must be greater than zero.")
#   }
# 
#   y <- rgamma(length(a), a, 1)
#   y / sum(y)
# }





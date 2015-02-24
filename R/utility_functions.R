# Title: utility_functions.R
# Author: Stuart Lee, Katherine Smith
# Description: Utility functions for the haldane package

#--- Stastical functions
#' Harmonic Mean
#'
#' This function returns the harmonic mean of an array.
#' Code from http://economistatlarge.com/r-guide/arithmetic-harmonic-geometric-means-r
#' @param array numeric vector
#' @keywords harmonic mean
#' @export
#' @examples
#' harmonicMean(rep(100,10))
harmonicMean <- function(array) {


  if (!is.numeric(array)) {
    stop("Passed argument must be an array. Consider using sapply for data frames.")
  }

  if (any(array<=0)) {
    stop("All values must be greater than zero.")
  }
  length(array) / sum(1 / array)
}

dirichletAlpha <- function(K,m,a0) {
  #' Concentration parameters for Dirichlet distribution
  #'
  #' @description  Generate vector alpha giving the parameters of a Dirichlet distribution
  #' @details A starting value a0 is multiplied by a factor m to obtain K categories
  #' @param K number of categories
  #' @param m weighting factor
  #' @param a0 initial value
  #' @export
  #'
  a0*m^(0:(K-1))
}


#' Dirichlet distribution random number generator
#'
#' @description Draw random numbers from a Dirichlet distribution with parameters given by the vector a.
#' Code modified from http://stats.stackexchange.com/a/30972
#'
#' @param a vector of concentration parameters
#' @return length(a) dirichlet random numbers
#' @keywords Dirichlet distribution
#' @export
#' @examples
#' rdirichlet(c(1,2,4))
rdirichlet <- function(a) {


  if (!is.numeric(a)) {
    stop("Argument must be a numeric vector.")
  }

  if (any(a<=0)) {
    stop("All values must be greater than zero.")
  }

  y <- rgamma(length(a), a, 1)
  y / sum(y)
}


# n is the number of SNPs S, m is the rate mut_rate, and t is the maximum value, to be set to 1
# rtexp(n=S, m=mut_rate, t=1)
# this code is modified from that given at https://stat.ethz.ch/pipermail/r-help/2010-June/242615.html
# look into alternate distrubition functions
#' Truncated exponential distribution random number generation
#'
#' @description Simulate random numbers from a truncated exponential distribution
#' @details Simulate SNP alternate allele frequencies from a truncated exponential distribution
#' Code is modified from that given at https://stat.ethz.ch/pipermail/r-help/2010-June/242615.html
#' @param n number of random values to return
#' @param mu rate parameter for exponential
#' @param t cut-off threshold for truncation
#' @return n random values from TEXP
#' @export
rtexp <- function(n, mu, t) {



  # first compute the inverse CDF for the truncated exponential. u is the quantile
  itexp <- function(u, mu, t) {
    -log(1 - u*(1 - exp(-t * mu)))/mu
  }
  x <- itexp(runif(n), mu, t)
  # assert that values lie in [0,t]
  stopifnot( min(x) >= 0, max(x) <= t )
  x
}

#--- Input/Output functions

#' Convert VCF file into Matrix format
#'
#' @description Read in Variant Call format file and output a matrix of SNPs
#' @details Assumes VCF file only contains SNPs and coming from haploid organism
#' . The output is n by m matrix, where n corresponds to sample name and
#' m is the number of SNPs. Each element of the matrix can take values 0, 1, 2
#' or 3 mapping to the set of genotypes {./., A/A, B/B, A/B} where A is reference
#' allele and B is alternate allele.
#' @param vcf vcf file to be read
#' @return a matrix of SNP calls
#' @export
vcftoMatrix <- function(vcf) {
  # establish connection to file
  # check whether gzip or not
  if (!grepl(".vcf", vcf)) stop("Not a valid vcf file.")
  if (grepl("vcf.gz$", vcf)) {
    con <- gzfile(vcf)
  }
  else {
    con <- file(vcf)
  }

  headerLine <- "^#"
  variants <- readLines(con, n = 1000)
  noheader <- grep(headerLine, variants)
  variants[noheader]

}
}

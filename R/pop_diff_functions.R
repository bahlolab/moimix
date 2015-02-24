# Title: pop_diff_functions.R
# Author: Stuart Lee, Katherine Smith
# Description: Functions for computing population genetic metrics for
# differentiating between populations.
# At the moment we have:
#   pi diversity
#   Jost's D
#   Nei's Gst
# TODO:
#   Fst
#
#
#

# calculate nucleotide diversity from an N x S matrix of alleles
piDiversity <- function (x) {
  #' Calculate the SNP diversity within a population
  #'
  #' @description Nucleotide diversity from a matrix of alleles
  #' @param x A matrix of size n.samples by n.snps
  #' @return numeric estimate of diversity
  #' @export
  # calculate the number of differences between all pairwise combinations of rows

  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)

  pairwiseDiff <- function(y) {
    rowSums(abs(x != matrix(y, ncol = nc, nrow = nr, byrow = TRUE)))
  }

  res <- apply(x, 1, pairwiseDiff)
  # the result is a symmetric matrix, so set lower triangle to NA to avoid double counting
  res[lower.tri(res)] <- NA

  # sum the number of nucleotide differences and weight by the frequency of each sequence and the number of variants
  nt.diffs <- sum(res, na.rm=T)
  pi.est <- nt.diffs/(nc * nr^2)
  pi.est
}

# calculate Jost's D
D_Jost <- function(n, H_T_est, H_S_est) {
  num <- n*(H_T_est - H_S_est)
  denom <- (n - 1)*(1 - H_S_est)
  return( num/denom )
}

# Calculate Nei's $G_{ST}$
G_ST <- function(H_T_est, H_S_est) {
  return( (H_T_est - H_S_est) / H_T_est )
}

# Calculate Jost's D and $G_{ST}$ using the equations on p4022 of Jost:2008 (DOI 10.1111/j.1365-294X.2008.03887.x)
# $2 \tilde(N)$ becomes $\tilde(N)$ because we have haploid data
# the input is a list of matrices of dimension C x S -- one for each deme
partition.diversity<-function(x) {

  # estimated allele frequencies
  # this is a matrix of dimension S x n.demes
  AAFs.est<-sapply(x, function(x) { colSums(x)/nrow(x) } )

  # calculate the 'heterozygosity' of each deme
  H_j <- 1 - AAFs.est^2 - ( 1 - AAFs.est)^2

  # H_S is the mean 'heterozygosity' of the individual demes
  H_S <- rowMeans( H_j )

  # the heterozygosity of the pooled subpopulations
  H_T <- 1 - (rowMeans(AAFs.est))^2 - (rowMeans(1 - AAFs.est))^2

  # harmonic mean of deme sample size
  N.harm <- harmonicMean( sapply(x,nrow) )

  # a nearly unbiased estimator of H_S
  ploidy<-1
  H_S_est <- ( H_S * ploidy * N.harm ) / ( ploidy * N.harm - 1 )

  # a nearly unbiased estimator of H_T
  H_T_est <- H_T + H_S_est/(ploidy * N.harm * n.demes)

  # calculate $G_{ST}$
  G.ST <- G_ST(H_T_est, H_S_est)

  # calculate $D_{Jost}$
  D.Jost <- D_Jost( length(x), H_T_est, H_S_est)

  return( list(
    N.harm = N.harm,
    AAFs.est = AAFs.est,
    H_j = H_j,
    H_S = H_S,
    H_T = H_T,
    G_ST = G.ST,
    D_Jost = D.Jost
  ))
}




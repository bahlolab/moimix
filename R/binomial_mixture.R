# binomial_mixture.R
# Author: Stuart Lee
# Date: 07/04/2015
# Description: Functions for clustering variant frequencies
# using a two component mixture of binomial distributions

#' Collapse allele frequency to half axis
#' @param p vector of proprotions
#' @export
collapseProbs <- function(p) { ifelse(p > 0.5, 1-p, p)}

#' return pdf of two component binomial mixture
dbinommix <- function(x, n1, n2, p, q, f) {
  m1 <- dbinom(x, size = n1, prob = p)
  m2 <- dbinom(x, size = n2, prob = q)
  f * m1 + (1-f) * m2
}
#' Compute complete log-likelihood based on current estimates
eStep <- function(x, theta) {

}

#' Maxisime likelihood given updated estimates
mStep <- function(theta.new) {

}

#' EM algorithm for two-component mixture model
#' @param x observed vector of data
#' @param epsilon precsion (default 0.01)
#' @param niter number of iterations (default 1000)
binommixEM <- function(x, epsilon = 0.01, niter = 1000) {

  # initialise parameters
  p.init <- runif(1)
  q.init <- runif(1)
  f.init <- runif(1)
  notdone <- TRUE
  while(notdone) {

    loglik <- eStep(x, theta)
    theta_hat <- mStep(loglik)

    error <- theta_hat - theta
    theta <- theta_hat
    if (all(abs(error) <= epsilon) && niter) {
      break
    }
    else {
      niter <- niter - 1
      if (niter == 0) {
        stop("EM algorithm did not converge.")
      }
    }
  }
}

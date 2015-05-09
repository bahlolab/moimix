# binomial_mixture.R
# Author: Stuart Lee
# Date: 07/04/2015
# Description: Functions for clustering variant frequencies
# using a two component mixture of binomial distributions

#' Collapse allele frequency to half axis
#' @param p vector of proprotions
#' @export
collapseProbs <- function(p) { ifelse(p > 0.5, 1-p, p)}

#' Update class probabilities (i.e compute tau)
updateClassProb <- function(x, k, mixture.comp, mixture.weights) {
  classProb <- function(x, i) {
    mixture.weights[i] * dbinom(x, size = 30, prob = mixture.comp[i])
  }
  probs <- sapply(1:k, FUN = function(k) classProb(x, k))
  probs / rowSums(probs)
}

updateWeights <- function(class.probs) {
  colMeans(class.probs)
}

updateComponents <- function(x, class.probs) {
  colSums(x * class.probs) / colSums(30 * class.probs)
}

mixLL <- function(x, k, class.memberships, mixture.weights, mixture.comp) {
  within.class.ll <- sapply(1:k,
         FUN = function(k) log(mixture.weights[k])*
           dbinom(x[class.memberships == k], size = 30, prob = mixture.comp[k], log = TRUE))
  sum(unlist(within.class.ll))
}

#' return pdf of two component binomial mixture
dbinommix <- function(x, n1, n2, p, q, f) {
  m1 <- dbinom(x, size = n1, prob = p)
  m2 <- dbinom(x, size = n2, prob = q)
  f * m1 + (1-f) * m2
}

#' EM algorithm for k-component binomial mixture model
#' @param x observed vector of data (read counts)
#' @param k number of mixture components (max=5)
#' @param epsilon precsion (default 1e-6)
#' @param mixture.comp vector of mixture paramters (length k)
#' @param mixture.weights vector of mixture weights (length k)
#' @param niter number of iterations (default 1000)
binommixEM <- function(x, k, mixture.comp = NULL, mixture.weights = NULL,
                       epsilon = 1e-6, niter = 1000) {

  # initialise parameters
  n = length(x)

  # if mixture model parameters not supplied
  # then guess them
  if (is.null(mixture.comp) && is.null(mixture.weights)) {
    initialiseEM()
  }

  # initialse log-likelihoods
  ll <- 0
  oldll <- 0
  notdone <- TRUE
  while(notdone) {
    oldll <- ll

    # update cluster memberships (eStep)
    class.probs <- updateClassProb(x, k, mixture.comp, mixture.weights)

    # mStep
    # update class memberships
    zhat <- apply(class.probs, 1, which.max)

    # update mixture proportions
    mixture.weights <- updateWeights(class.probs)

    # update binomial probabilites
    mixture.comp <- updateComponents(x, class.probs)

    # compute new logliklihood
    ll <- mixLL(x, k, zhat, mixture.weights, mixture.comp)
    error <- theta_hat - theta
    theta <- theta_hat
    if ((abs(ll - oldll) <= epsilon) && niter) {
      message(paste("EM algorithm converged at ", niter))
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

# binomial_mixture.R
# Author: Stuart Lee
# Date: 07/04/2015
# Description: Functions for clustering variant frequencies
# using a two component mixture of binomial distributions

#' Collapse allele frequency to half axis
#' @param p vector of proprotions
#' @export
collapseProbs <- function(p) { ifelse(p > 0.5, 1-p, p)}

#' Update class probabilities (i.e compute tau in Estep)
updateClassProb <- function(x, N, k, mixture.comp, mixture.weights) {
  classProb <- function(x, i) {
    mixture.weights[i] * dbinom(x, size = N, prob = mixture.comp[i])
  }
  probs <- sapply(1:k, FUN = function(k) classProb(x, k))
  probs / rowSums(probs)
}

#' update class weights
updateWeights <- function(class.probs) {
  colMeans(class.probs)
}
#' update probs for binomials
updateComponents <- function(x, N, class.probs) {
  colSums(x * class.probs) / colSums(N * class.probs)
}

#' Likelihood function for mixture model
mixLL <- function(x, N, k, class.memberships, mixture.weights, mixture.comp) {
  within.class.ll <- sapply(1:k,
         FUN = function(k) log(mixture.weights[k])*
           dbinom(x[class.memberships == k],
                  size = N[class.memberships == k],
                  prob = mixture.comp[k],
                  log = TRUE))
  sum(unlist(within.class.ll))
}

# #' return pdf of two component binomial mixture
# dbinommix <- function(x, N, p, q, f) {
#   m1 <- dbinom(x, size = n1, prob = p)
#   m2 <- dbinom(x, size = n2, prob = q)
#   f * m1 + (1-f) * m2
# }

#' Initialise mixture model parameters for EM
#' uses kmeans on proportions
initEM <- function(x, N, k) {
  # use kmeans with k clusters
  kfit <- kmeans(x/N, k, nstart = 20)
  mixture.comp <- kfit$centers
  mixture.weights <- as.numeric(table(kfit$cluster))/length(x)
  list(mixture.comp = mixture.comp,
       mixture.weights = mixture.weights)
}
#' EM algorithm for k-component binomial mixture model
#' @param x observed vector of counts (read counts supporting SNV)
#' @param N observed vector of counts (coverage at SNV)
#' @param k number of mixture components (max=5)
#' @param epsilon precsion (default 1e-6)
#' @param mixture.comp vector of mixture paramters (length k)
#' @param mixture.weights vector of mixture weights (length k)
#' @param niter number of iterations (default 1000)
#' @export
binommixEM <- function(x, N, k, mixture.comp = NULL, mixture.weights = NULL,
                       epsilon = 1e-6, niter = 1000) {

  # input error checking
  # number of components out of bounds
  if (k < 1 | k > 5) {
    stop("Number of mixture components must be between 1 and 5")
  }
  # mixture model parameters missing
  if ((is.null(mixture.comp) & !is.null(mixture.weights)) |
        (!is.null(mixture.comp) & is.null(mixture.weights))) {
    stop("Mixture model parameters must both be NULL or both initialised")
  }
  # input parameters not same length as k
  if (!(is.null(c(mixture.comp, mixture.weights))) &&
        ((length(mixture.comp) < k) | (length(mixture.weights) < k))) {
    stop(paste("Mixture model parameters must have length:", k))
  }

  # reads and coverage vectors not same length
  if (length(x) != length(N)) {
    stop("X and N must have same length")
  }
  # initialise parameters
  n = length(x)

  # if mixture model parameters not supplied
  # then guess them using k-means
  if (is.null(c(mixture.comp,mixture.weights))) {
    init.params <- initEM(x, N, k)
    mixture.comp <- init.params$mixture.comp
    mixture.weights <- init.params$mixture.weights
  }

  # initialse log-likelihoods
  ll <- 0
  oldll <- 0

  print("Start EM algorithm")
  print(paste("The parameters (mu, pi): ",
              paste(c(mixture.comp, mixture.weights), collapse=", ")))

  nstart <- niter
  while(TRUE) {
    oldll <- ll

    # update cluster memberships (eStep)
    class.probs <- updateClassProb(x, N, k, mixture.comp, mixture.weights)

    # mStep
    # update class memberships
    zhat <- apply(class.probs, 1, which.max)

    # update mixture proportions
    mixture.weights <- updateWeights(class.probs)

    # update binomial probabilites
    mixture.comp <- updateComponents(x, N, class.probs)

    # compute new logliklihood
    ll <- mixLL(x, N, k, zhat, mixture.weights, mixture.comp)
    print(paste("The log-like is:", ll))
    if ((abs(ll - oldll) <= epsilon) && niter) {
      message(paste("EM algorithm converged after ",
                    nstart - niter, "iteration"))
      break
    }
    else {
      niter <- niter - 1
      if (niter == 0) {
        stop("EM algorithm did not converge.")
      }
    }
  }
  # update cluster memberships based on final update
  cluster.probs <- updateClassProb(x, N, k, mixture.comp, mixture.weights)
  cluster.memberships <- apply(cluster.probs, 1, which.max)
  return(list(cluster.memberships = cluster.memberships,
              cluster.probs = cluster.probs,
              pi = mixture.weights,
              mu = mixture.comp,
              log.lik = ll))

}

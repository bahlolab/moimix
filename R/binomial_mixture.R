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
  #classProb <- function(x, i) {
  #  mixture.weights[i] * dbinom(x, size = N, prob = mixture.comp[i])
  #}
  probs <- sapply(1:k, FUN = function(j) 
      mixture.weights[j] * dbinom(x, size = N, prob = mixture.comp[j]))
  denom <- rowSums(probs)
  tau <- probs / denom
  tau  
}

#' update class weights
updateWeights <- function(class.probs) {
    colMeans(class.probs)
}
#' update probs for binomials
updateComponents <- function(x, N, class.probs) {
  colSums(class.probs*x) / colSums(class.probs*N)
}

#' log-Likelihood function for mixture model
# mixLL <- function(x, N, k, mixture.comp, mixture.weights) {
#     # generate log(pi_k) + log(dbinom(x, N, mu_k))
#     # gives us a length(x) by k
#     within.class.ll <- sapply(1:k,
#                               function(i) {
#                                   log(mixture.weights[i]) +
#                                       dbinom(x, size = N, prob=mixture.comp[i], log = TRUE)
#                               } )
#     print(mixture.comp)
#     print(mixture.weights)
#     #print(within.class.ll)
#     #print(class.probs)
#     # we also have a length(x) by k responsibilites matrix
#     # take dot product of rows with columns then sum to get
#     # expected log-likelihood
#     sum(rowSums(within.class.ll*class.probs))
#     
# }
#' Expected log-Likelihood function for mixture model
mixLL <- function(x, N, k, class.probs, mixture.comp, mixture.weights) {
    # generate log(pi_k) + log(dbinom(x, N, mu_k))
    # gives us a k by length(x) matrix
    within.class.ll <- apply(sapply(mixture.comp, dbinom,
                                    x = x,
                                    size = N,
                                    log = TRUE), 1, 
                             function(i) i + log(mixture.weights))
    # we also have a length(x) by k responsibilites matrix
    # take dot product of rows with columns then sum to get
    # expected log-likelihood
    sum(sapply(1:length(x), 
               function(i) sum(class.probs[i,] * within.class.ll[,i])))
    
}


#' Initialise mixture model parameters for EM
#' uses kmeans on proportions
initEM <- function(x, N, k) {
  # use kmeans with k clusters
  kfit <- kmeans(x/N, k, nstart = 20)
  mixture.comp <- as.numeric(kfit$centers)
  mixture.weights <- tabulate(kfit$cluster)/length(x)
  list(mixture.comp = sort(mixture.comp),
       mixture.weights = sort(mixture.weights))
}
#' EM algorithm for k-component binomial mixture model
#' @param x observed vector of counts (read counts supporting SNV)
#' @param N observed vector of counts (coverage at SNV)
#' @param k number of mixture components (max=5)
#' @param epsilon precsion (default 1e-6)
#' @param mixture.comp vector of mixture paramters (length k)
#' @param mixture.weights vector of mixture weights (length k)
#' @param niter number of iterations (default 1000)
#' @param verbose print running of EM algorithm (default FALSE)
#' @return a list with the following objects
#'  n total number of SNVs
#'  k number of mixture components
#'  cluster.memberships assignments to each class
#'  cluster.probs an n by k matrix with probability of SNV membership to each class
#'  pi a k length vector with estimated mixture weights
#'  mu a k length vector with estimated mixture components
#' @export
binommixEM <- function(x, N, k, mixture.comp = NULL, mixture.weights = NULL,
                       epsilon = 1e-6, niter = 1000, verbose = FALSE) {

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
    if (length(N) == 1) {
      # range checking
      stopifnot(0 <= min(x) && max(x) <= N)
      if(verbose) {
        message(paste0("Assuming uniform coverage: ", N))
      }

    }
    else {
      stop("X and N must have same length")
    }
  }
  stopifnot(0 <= min(x) && max(x) <= max(N))

  # initialise parameters
  n = length(x)

  # if mixture model parameters not supplied
  # then guess them using k-means
  if (is.null(c(mixture.comp,mixture.weights))) {
      init.params <- initEM(x, N, k)
      mixture.comp <- init.params$mixture.comp
      mixture.weights <- init.params$mixture.weights
  }
  else {
      # sort paramters for identifiability
      mixture.comp <- sort(mixture.comp)
      mixture.weights <- sort(mixture.weights)
  }

  # initialse log-likelihoods
  ll <- 0
  oldll <- 0

  if(verbose) {
    print("Start EM algorithm")
    print("The initial parameter guesses are (pi, mu):")
    print(c(mixture.weights, mixture.comp))
  }

  nstart <- niter
  while(TRUE) {
    
    oldll <- ll
      
    # update cluster memberships (eStep)
    class.probs <- updateClassProb(x, N, k, mixture.comp, mixture.weights)

    # mStep
    # update mixture proportions
    mixture.weights <- updateWeights(class.probs)
    # update binomial probabilites
    mixture.comp <- updateComponents(x, N, class.probs)
    # recompute expected log-likelihood with current guesses
    ll <- mixLL(x, N, k, class.probs, mixture.comp, mixture.weights)

    if (verbose) { print(paste("The log-like is:", ll)) }

    if ((abs(ll - oldll) <= epsilon) && niter) {
      if(verbose){
        message(paste("EM algorithm converged after ",
                      nstart - niter, "iteration"))
      }
      break
    }
    else {
      niter <- niter - 1
      if (niter == 0) {
        message("EM algorithm did not converge.")
        break
      }
    }

  }
  convergence.iter <- (nstart - niter)
  convergence.alg <- (nstart > 0)
  # update cluster memberships based on final update
  cluster.probs <- updateClassProb(x, N, k, mixture.comp, mixture.weights)
  cluster.memberships <- apply(cluster.probs, 1, which.max)
  return(list(n = n,
              k = k,
              cluster.memberships = cluster.memberships,
              cluster.probs = cluster.probs,
              pi = mixture.weights,
              mu = mixture.comp,
              log.lik = ll,
              converge.iter = convergence.iter,
              converge.true = convergence.alg))

}

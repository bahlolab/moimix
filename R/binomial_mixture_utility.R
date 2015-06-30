# Title: binomial_mixture_utility.R
# Description: Tools for model selection and evalution
# for binomial mixture model.
# Author: Stuart Lee
# Date: 11/05/2015

#' Bayesian Information Criterion for binomial mixture model
#'
#' @param mixture estimated mixture model
#' @return The BIC of a model
#' @export
bic <- function(x, N, mixture) {
  # we multiply k by 2 since we estimate
  # 2k paramters in the binommix
  2*mixture$k*log(mixture$n) -2 * mixLL(x, N, mixture$k, mixture$pi, mixture$mu)

}

#' Akaike Infomation Criterion for binomial mixture model
#'
#' @inheritParams bic
#' @return The AIC of a model
#' @export
aic <- function(x, N, mixture) {
  2 * (2*mixture$k - mixLL(x, N, mixture$k, mixture$pi, mixture$mu))
}

#' Mean Square Error for binomial mixture model
#'
#' @param mixture estimated mixture model
#' @param theta     true mixture-model parameters of form (pi_1,..,pi_k, mu_1,...mu_k)
#' @export
mseMM <- function(mixture, theta) {
  # potential matching problem
  # which we avoid by sorting the components
  stopifnot(length(theta) == 2*mixture$k)
  stopifnot(sum(theta[1:mixture$k]) == 1)

  k.index <- 1:mixture$k
  pi.hat <- sort(mixture$pi)
  mu.hat <- sort(mixture$mu)
  pi <- sort(theta[1:k])
  mu <- sort(theta[(k+1):(2*k)])
  print(c(pi, pi.hat))
  print(c(mu, mu.hat))
  mse <- list(pi.mse = mean((pi - pi.hat)^2),
              mu.mse = mean((mu - mu.hat)^2),
              all.mse = mean((c(pi,mu) - c(pi.hat,mu.hat))^2))
  return(mse)
}

#' Use an information criterion to perform model selection
#' for number of mixture components
#'
#' @param x vector of observed read counts
#' @param N coverage
#' @param k maximum number of components to estimate
#' @param method string with information criterion either 'bic' or 'aic'
#' @param ... other parameters passed to binommixEM
#' @export
autoSelect <- function(x, N, k, method, ...) {
  if(k < 1 | k > 5) {
    stop("Maximum number of components must be between 1 and 5")
  }
  if(!(method %in% c("aic", "bic"))) {
    stop("Method must either be aic or bic")
  }
  if(method == "aic") {
    selector <- haldane::aic
  }
  if(method == "bic") {
    selector <- haldane::bic
  }

  # initalise
  crits <- c()
  models <- list()
  for (i in 1:k) {
    model <- binommixEM(x, N, i, ...)
    ic <- selector(x, N, model)
    models[[i]] <- list(mu = sort(model$mu), pi=sort(model$pi))
    crits[i] <- ic
  }

  min.ic <- which.min(crits)

  results <- data.frame(component = 1:k, values = crits)
  best.model <- models[[min.ic]]

  return(list(info.crit = results,
              k = min.ic,
              mu = best.model$mu,
              pi = best.model$pi))
}

#' CDF for binomial mixture model
#'
#'@param mixture
#'@param x vector of read counts supporting each SNV
#'@param N vector of coverage at SNV
#'@return Vector of theoretical quantiles
binommixCDF <- function(mixture, x, N) {
  pi <- mixture$pi
  k <- mixture$k
  pbinomForMix <- function(x, N, component) {
    pi[component] * pbinom(x, size = N,
                           prob = mixture$mu[component])
  }

  pbinoms <- sapply(1:k, pbinomForMix, x = x, N = N)
  return(rowSums(pbinoms))
}

#' Fisher Information Approximation
#'
#' @description Estimate the Fisher information matrix for
#' a bionimal mixture model
#'
#' @param x vector of read counts
#' @param N vector of coverage at site
#' @param k number of components
#' @param mu mixture components
#' @param pi mixture weights
#' @return A 2k by 2k Fisher information matrix
#' @export
infoMat <- function(x, N, k, mu, pi) {

  class.probs <- updateClassProb(x, N, k, mu, pi)
  mat <- matrix(0, ncol = 2*k, nrow = 2*k)
  # -- second-derivative evaluated at mixture weights
  dq2dpi2 <- function(class.probs, x, N, pi, i) {
    -sum(class.probs[,i]/(pi[i]^2)) - sum((1 - class.probs[,i]) / (1 - pi[i])^2)
  }
  # -- second-derivative evaluated at mixture components
  dq2dmu2 <- function(class.probs, x, N, mu, i) {
    -sum((class.probs[, i] * x )/ (mu[i])^2) -sum(class.probs[,i]*(N-x)/(mu[i])^2)
  }
  diag(mat)[1:k] <- sapply(1:k,
                           function(i)
                             dq2dpi2(class.probs, x, N, pi, i))

  diag(mat)[(k+1):(2*k)] <- sapply(1:k,
                                   function(i)
                                     dq2dmu2(class.probs, x, N, mu, i))

  return(-mat)
}



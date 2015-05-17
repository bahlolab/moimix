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
bic <- function(mixture) {
  # we multiply k by 2 since we estimate
  # 2k paramters in the binommix
  -2 * mixture$log.lik + 2*mixture$k*log(mixture$n)
}

#' Akaike Infomation Criterion for binomial mixture model
#'
#' @inheritParams bic
#' @return The AIC of a model
#' @export
aic <- function(mixture) {
  2 * (2*mixture$k - mixture$log.lik)
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
#' @export
autoSelect <- function(x, N, k, method) {
  if(k < 1 | k > 5) {
    stop("Maximum number of components must be between 1 and 5")
  }
  if(!(method %in% c("aic", "bic"))) {
    stop("Method must either be aic or bic")
  }
  selector <- match.fun(method)

  # initalise
  crits <- c()
  models <- list()
  for (i in 1:k) {
    model <- binommixEM(x, N, i)
    ic <- selector(model)
    models[[i]] <- list(model = model)
    crits[i] <- ic
  }

  min.ic <- which.min(crits)

  results <- data.frame(component = 1:k, values = crits)
  best.model <- model[[min.ic]]

  return(list(method = method,
              info.crit = results,
              model = best.model))
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



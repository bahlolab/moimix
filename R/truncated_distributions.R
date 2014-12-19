# Title: truncated_distributions.R
# Description: Code to generate truncated distrubtion functions
# for a beta distribution and poisson distribution
# see paper: http://www.jstatsoft.org/v16/c02/paper for details
# Author: Stuart Lee
# 30/09/2014

# --- Beta functions
#' Truncated beta distribution
#'
#' @description Density, distribution function, quantile function and
#' random generation for the truncated Beta distribution with parameters
#' shape1 and shape2 (and optional truncation bounds)
#'
#' @param x,q vector of quantiles
#' @param shape1,shape2 non-negative parameters of Beta distributions
#' @param lower,upper bounds for truncation
#' @export
dtbeta <- function(x, shape1, shape2, lower = 0, upper = 1) {

  stopifnot(lower <= upper,shape2 > 0, shape1 > 0)
  tt <- rep(0, length(x))
  normalize.factor <- pbeta(upper, shape1, shape2) - pbeta(lower, shape1, shape2)
  tt[x >= lower & x <= upper] <- dbeta(x[x >= lower & x <= upper],
                                       shape1, shape2) / normalize.factor
  return(tt)
}

#' @describeIn dtbeta Cumulative distribution function
#' for truncated beta
ptbeta <- function(x, shape1, shape2, lower = 0, upper = 1) {
  stopifnot(lower <= upper,shape2 > 0, shape1 > 0)
  tt <- x
  aa <- rep(lower, length(x))
  bb <- rep(upper, length(x))
  normalize.factor <- pbeta(bb, shape1, shape2) - pbeta(aa, shape1, shape2)
  tt <- pbeta(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), shape1, shape2)
  tt <- (tt - pbeta(aa, shape2, shape1)) / normalize.factor
  return(tt)
}

#'  @describeIn dtbeta Quantile function
#'  @param p vector of probabilities
qtbeta <- function(p, shape1, shape2, lower = 0, upper = 1) {
  tt <- p
  pin <- pbeta(lower, shape1, shape2) + p * (pbeta(upper, shape1, shape2) - pbeta(lower, shape1, shape2))
  tt <- qbeta(pin, shape1, shape2)
  return(tt)
}

#' @describeIn dtbeta  Random number generation
#' @param n number of observations
rtbeta <- function(n, shape1, shape2, lower = 0, upper = 1) {
  u <- runif(n, min = 0, max = 1)
  x <- qtbeta(u, shape1, shape2, lower, upper)
  return(x)
}

# -- likelihood and methods of moments function
beta.mom <- function(x, lower = 0.01, upper = 100) {
  # method of moments for a Beta distribution
  x.bar <- mean(x)
  n <- length(x)
  v <- var(x) * (n - 1)/n
  R <- 1/x.bar - 1

  f <- function(a) {
    # note: undefined when a=0
    R * a^2/((a/x.bar)^2 * (a/x.bar + 1)) - v
  }

  u <- uniroot(f, c(lower, upper))

  return(c(shape1 = u$root, shape2 = u$root * R))
}

# log-likelihood for truncated distribution
lltrunc_beta <- function(shape.par, samples, low = 0, high = 1) {
  shape1 = shape.par[1]
  shape2 = shape.par[2]
  likelihood <- dtrunc_beta(samples, shape1, shape2, low, high)
  -sum(log(likelihood))
}

dtpois <- function(xs, threshold, lambda) {
  ifelse(xs > threshold, dpois(xs, lambda) / ppois(threshold, lambda, lower.tail = FALSE), 0) }

# negative log-likelihood for truncated poisson
ll_dtpois <- function(xs, threshold, lambda) {
  -sum(log(dtpois(xs, threshold, lambda)))
}

# Title: binomial_mixture_utility.R
# Description: Tools for model selection and evalution
# for binomial mixture model.
# Author: Stuart Lee
# Date: 11/05/2015
#' Collapse allele frequency to half axis
#' @param p vector of proprotions
#' @export
collapseProbs <- function(p) { ifelse(p > 0.5, 1-p, p)}

#' Filter read-count matrix
#' @param x vector of b-allele depths
#' @param N vector of site depths
#' @param p vector containing boundary cut-offs for inclusion
#' @return matrix with two columns
#' @export
filterCounts <- function(x, N) {
    #if (length(p) != 2) stop("Boundaries proportion must be between 0 and 1")
    #if (sum(p) != 1) stop("p must sum to 1")
    # if (diff(p) == 0) stop("no data left!")
    diff <- N - x
    index = which(x  > 0)
    cbind(x[index], diff[index])
}

#' Convert odds to probability
toProb <- function(o) {
  exp(o) / (1+exp(o))
}

#' Mean Square Error for binomial mixture model
#'
#' @description Compute sum of squares for mixture components
#' @param mixture estimated mixture model
#' @param theta   true mixture-model parameters of form (pi_1,..,pi_k, mu_1,...mu_k)
#' @importMethodsFrom flexmix parameters prior
#' @export
mse <- function(model, theta) {
  stopifnot(is(model, "flexmix"))
  stopifnot(length(theta) == 2*model@k)

  k <- model@k
  estimates <- getTheta(model)
  estimates <- estimates[order(-estimates$pi.hat),]
  pi <- theta[1:k]
  # check pi's sum to 1
  stopifnot(all.equal(sum(pi), 1))
  
  mu <- theta[(k+1):(2*k)]
  true <- data.frame(pi = pi, mu = mu)
  # order by pi component
  true <- true[order(-true$pi),]
  print(estimates)
  print(true)
  # evaluate the error in each component
  mse <- data.frame(pi.mse = sum((true$pi - estimates$pi.hat)^2)/k,
              mu.mse = sum((true$mu - estimates$mu.hat)^2)/k)
  return(mse)
}

#' CDF for binomial mixture model
#'
#'@param model flexmix object
#'@param x vector of read counts supporting each SNV
#'@param N vector of coverage at SNV
#'@return Vector of theoretical quantiles
binommixCDF <- function(model, x, N) {
    estimates <- getTheta(model)
    pi <- estimates$pi.hat
    mu <- estimates$mu.hat
    k <- model@k
    pbinomForMix <- function(x, N, component) {
        pi[component] * pbinom(x, 
                               size = N,
                               prob = mu[component])
    }
    
    pbinoms <- sapply(1:k, pbinomForMix, x = x, N = N)
    return(rowSums(pbinoms))
}

#' Fisher Information Approximation
#'
#' @description Estimate the Fisher information matrix for
#' a bionimal mixture model
#'
#' @param y a two column matrix of read counts
#' @param model flexmix object
#' @return A 2k by 2k Fisher information matrix
#' @export
infomat <- function(y, model) {
    
    k <- model@k
    mu <- getTheta(model)$mu.hat
    pi <- getTheta(model)$pi.hat
    
    class.probs <- posterior(model)
    mat <- matrix(0, ncol = 2*k, nrow = 2*k)
    # -- second-derivative evaluated at mixture weights
    dq2dpi2 <- function(class.probs, pi, i) {
        -sum(class.probs[,i]) / (pi[i])^2
    }
    # -- second-derivative evaluated at mixture components
    dq2dmu2 <- function(class.probs, y, mu, i) {
        sum(class.probs[,i] * (-y[,1]/(mu[i]^2) + y[,2]/(1-(mu[i])^2)))
    }
    diag(mat)[1:k] <- sapply(1:k,
                             function(i)
                                 dq2dpi2(class.probs, pi, i))
    
    diag(mat)[(k+1):(2*k)] <- sapply(1:k,
                                     function(i)
                                         dq2dmu2(class.probs, y, mu, i))
    
    return(-mat)
}

#' Compute standard errors for mixture model parameter estimates
#' 
#' @description se is computed as the inverse of the Fisher information matrix
#' at the EM estimates
#' @param y a two column matrix of read counts
#' @param model flexmix object
#' @export
seMM <- function(y, model) {
    k <- model@k
    info.mat.est <- infomat(y, model)
    # compute inverse
    inv.info.mat.est <- solve(info.mat.est)
    # return diagonals
    se <- sqrt(diag(inv.info.mat.est))
    se.pi <- se[1:k]
    se.mu <- se[(k+1):(2*k)]
    data.frame(se.pi = se.pi, 
               se.mu = se.mu)
}

#' Generate Confidence Intervals for
#' parameter estimates using estimated Fisher Information
#' Matrix
#' 
#' @param y filtered counts matrix
#' @param model flexmix object
#' @param alpha confidence level (default = 0.05)
#' @return a data.frame containing parameter estimates
#' and lower and upper confidence limits
ciMM <- function(y, model, alpha = 0.05) {
    z <- c(qnorm(alpha/2), qnorm(1-alpha/2))
    bounds.est <- seMM(y, model)
    param.est <- getTheta(model)
    ci.est <- data.frame(pi.hat = param.est$pi.hat,
                         pi.lower = z[1] * bounds.est$se.pi,
                         pi.upper = z[2] * bounds.est$se.pi,
                         mu.hat = param.est$mu.hat,
                         mu.lower = z[1] * bounds.est$se.mu,
                         mu.upper = z[2] * bounds.est$se.mu)
    ci.est
}

#' Generate iid mixture random variables 
#' for testing EM implementation
#' 
#' @param n number of realisations
#' @param N number of trials
#' @param k number of components
#' @param mu true mixture components
#' @param pi true mixture weights
#' @export
sampleMM <- function(n, N, k, mu, pi) {
    # assert that pi and mu must have length k
    stopifnot(length(pi) == k)
    stopifnot(length(mu) == k)
    # generate true hidden states
    states <- sample.int(k, size = n, replace = TRUE, prob = pi)
    # sample from binomial according to states
    observations <- rbinom(n, size = N, prob = mu[states])
    list(obs = observations, states = states)
}

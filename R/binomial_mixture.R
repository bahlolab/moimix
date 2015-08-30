# binomial_mixture.R
# Author: Stuart Lee
# Date: 07/04/2015
# Description: Functions for clustering variant frequencies
# using a two component mixture of binomial distributions
#' Update class probabilities (i.e compute tau in Estep)
#' @export
updateClassProb <- function(x, N, k, mixture.comp, mixture.weights) {
    #classProb <- function(x, i) {
    #  mixture.weights[i] * dbinom(x, size = N, prob = mixture.comp[i])
    #}
    # asster that mixture weights must sum to 1
    stopifnot(all.equal(sum(mixture.weights), 1))
    # shift everything into log space to avoid underflow    
    probs <- sapply(1:k, FUN = function(j) 
        log(mixture.weights[j])  + dbinom(x, 
                                    size = N, 
                                    prob = mixture.comp[j],
                                    log = TRUE))
    # probs are in log space
    denom <- logSum(probs)
    tau <- probs - denom
    # assert that all row probs must sum to 1
    stopifnot(all.equal(rowSums(exp(tau)), rep(1, nrow(tau))))
    
    list(log.probs = probs, tau = tau)
}
#' Take log of sum for logs
#' @export
logSum <- function(probs) {
    # compute log(sum_i^k pi_k*b(ni, mu_k)) from
    # log(pi_k*b(ni, mu_k)))
    # find maxP for each row of the matrix
    max.index <- apply(probs, 1, which.max)
    
    sapply(1:nrow(probs),
           function(j) 
               probs[j, max.index[j]] + 
               log1p(sum(exp(probs[j, -max.index[j]] - probs[j, max.index[j]]))))

}

#' update class weights
#' @export
updateWeights <- function(tau) {
    colSums(tau) / sum(colSums(tau))
}
#' update probs for binomials
#' @export
updateComponents <- function(x, N, tau) {
    colSums(tau*x) / colSums(tau*N)
}

#' full mstep combine pi and mu
#' @export
mStep <- function(x, N, class.probs) {
    # take inverse (get out of log-space)
    tau <- exp(class.probs)
    mixture.comp <- updateComponents(x, N, tau)
    mixture.weights <- updateWeights(tau)
    list(mixture.comp = mixture.comp, 
         mixture.weights = mixture.weights)
    
}

#' Expected log-Likelihood function for mixture model
#' @export
mixLL <- function(x, N, k, class.probs, mixture.comp, mixture.weights) {
    # generate log(pi_k) + log(dbinom(x, N, mu_k))
    # assert that mixture weights sum to 1
    stopifnot(all.equal(sum(mixture.weights), 1))
    
    pi.mat <- matrix(rep(log(mixture.weights), length(x)), 
                     ncol = k, 
                     byrow = TRUE)
    # genereate within cluster log-likelihoods
    within.class.ll <- sapply(mixture.comp, dbinom,
                                    x = x,
                                    size = N,
                                    log = TRUE)
    # now we have the log terms together in length(x) by k
    # matrix
    tau <- exp(class.probs)
    
    all.ll <- tau * (pi.mat + within.class.ll)
    
    sum(all.ll)
    # we also have a length(x) by k responsibilites matrix
    # multiply rows then take sums to get
    # expected log-likelihood
    
    # sum(rowSums(tau * all.ll))
    
}

# mixLL <- function(log.probs, tau) {
#     sum(log.probs * exp(tau))
# }

#' Seeding methods for mixture model
#' kmeans initilisation
#' @export
kmeans_seed <- function(x, N, k) {
    p <- as.matrix(x/N)    
    kfit <- kmeans(p, k, nstart = 20)
    mixture.comp <- as.numeric(kfit$centers)
    mixture.weights <- tabulate(kfit$cluster)/length(x)
    list(mixture.comp = mixture.comp,
         mixture.weights = mixture.weights)
}

#' random start seeding
#' choose array of starting points choose iteration that
#' gives best estimate
#' @export
random_seed <- function(x, N, k, nstart = 20) {
    # set up ll matrix
    ll.max <- rep(NA, nstart)
    # if k = 1 all pi's are 1, just guess mu as MLE
    if (k == 1) {
        
        return(list(mixture.comp = sum(x) / (N * length(x)),
                    mixture.weights = 1))
    }
    else {
        # sample mixture weights from Dirichlet
        pi.guess <-  t(MCMCpack::rdirichlet(nstart, rep(1, k)))
        # randomly partition data according to mixture weights
        # estimate mean within each group
        kk.partition <- apply(pi.guess, 2,
                              function(p) 
                                  sample.int(k, 
                                             size = length(x), 
                                             replace = TRUE, 
                                             prob = p))
        
        mu.guess <- apply(kk.partition, 2,
                           function(sub) 
                               sapply(1:k,
                                      function(j)
                                          mean(x[sub==j]) / N))
    }
    
    # compute 1 iteration for each starting setting
    for(i in 1:nstart) {
        estep <- updateClassProb(x, N, k, 
                                 mixture.comp = mu.guess[,i],
                                 mixture.weights = pi.guess[,i])        
        ll.max[i] <- mixLL(x, N, k, estep$tau, mu.guess[,i], pi.guess[,i])
    }
    
    # best start
    max.index <- which.max(ll.max)
    # sort for identifiability
    pi.best <- sort(pi.guess[, max.index])
    mu.best <- mu.guess[order(pi.guess[, max.index]), max.index] 
    list(mixture.comp = mu.best,
         mixture.weights = pi.best)
    
}

# grid search start
grid_seed <- function(x, N, k) {
    
}

#' Initialise mixture model parameters for EM
#' uses random start by default
#' @export
initEM <- function(x, N, k, method = "random_seed") {
    seed.fun <- match.fun(method)
    seed.fun(x, N, k)
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

  # assert mixture weiths sum to 1
  stopifnot(!is.null(mixture.weights) && sum(mixture.weights) == 1)
  
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
      # assert that mixture weights sum to 1
      stopifnot(sum(mixture.weights) == 1)
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
    
    current <- list(mixture.comp = mixture.comp,
                    mixture.weights = mixture.weights)
    # update cluster responsibilities (eStep)
    # print out is logged to avoid underflow
    estep <- updateClassProb(x, N, k, 
                             current$mixture.comp, 
                             current$mixture.weights)

 
    # mStep - new guesses
    theta <- mStep(x, N, estep$tau)
    mixture.weights <- theta$mixture.weights
    mixture.comp <- theta$mixture.comp
    
    

    # recompute expected log-likelihood with new guesses
    ll <- mixLL(x, N, k, estep$tau, mixture.comp, mixture.weights)    
    
    
    if (verbose) {
        print(paste("Current estimates are ", 
                    paste(c(mixture.weights, mixture.comp), collapse = ", "))) 
        print(paste("The log-like is:", ll)) 
    }

    # assert thtat log-likelihood increasing
    if((ll < oldll) && (niter < nstart)) {
        print(c(ll, oldll))
        stop("Log-likelihood is decreasing, but should increase at
             every iteration")
    }
    
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
  convergence.alg <- (niter != 0)
  pi <- mixture.weights
  mu <- mixture.comp
  
  return(list(n = n,
              k = k,
              pi = pi,
              mu = mu,
              log.lik = ll,
              converge.iter = convergence.iter,
              converge.true = convergence.alg))

}

binommix <- function(x, N, k, niter = 1000) {
  require(flexmix)
  y <- cbind(x, N-x)
  initFlexmix(y ~ 1, 
              k = k, 
              model = FLXMRglm(y ~ ., family = "binomial"),
              control = list(iter.max = niter,
                             minprior = 0))
}

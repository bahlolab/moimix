# binomial_mixture.R
# Author: Stuart Lee
# Date: 07/04/2015


#' Fit binomial mixture model on coverage data
#' 
#' @importFrom flexmix initFlexmix FLXMRglm
#' @param y matrix with 2 columns, first column contains 
#' counts in support of SNV, second in support of ref
#' @param N vector of depth at SNV sites
#' @param k vector of mixture components to fit
#' @param niter number of iterations to run 
binommix <- function(y, k, niter = 1000) {
    # I/O error handling
    if(!is(y, "matrix")) stop("y must be a matrix of counts")
    if(dim(y)[2] != 2) stop("y must have two columns")
    if (any(k < 1 | k > 5)) stop("Number of mixture components must be between 1 and 5")
    
    flexmix::initFlexmix(y ~ 1, 
                         k = k, 
                         model = flexmix::FLXMRglm(y ~ ., family = "binomial"),
                         control = list(iter.max = niter,
                                        minprior = 0))
}

#' Return estimated model parameters
#' 
#' @importFrom flexmix getModel parameters prior
#' @param model stepFlexmix or flexmix object
#' @param k if stepFlexmix choose model with k components
#' @param criterion if stepFlexmix choose model according to information 
#' criterion
getTheta <- function(model, k = NULL, criterion = NULL) {
    # error handling
    if (is(model, "stepFlexmix")) {
        if(is.null(k) && is.null(criterion)) {
            stop("model is stepFlexmix object, 
                 provide information criterion or number of components")
        }
        
        # model selection for stepMix 
        if(!is.null(criterion)) {
            stopifnot(criterion %in% c("AIC","BIC", "ICL"))
            model.select <- getModel(model, criterion)
        }
        else if (!is.null(k)) {
            model.select <- getModel(model, k)
        }
        data.frame(pi.hat = prior(model.select),
                   mu.hat = toProb(parameters(model.select)),
                   row.names = NULL)
    }
    else if(is(model, "flexmix")) {
        data.frame(pi.hat = prior(model),
                   mu.hat = toProb(parameters(model)),
                   row.names = NULL)
    }
    else {
        stop("model must be flexmix or stepFlexmix class")
    }
}

# binomial_mixture.R
# Author: Stuart Lee
# Date: 07/04/2015

#' Fit binomial mixture model on coverage data
#' 
#' @importFrom flexmix initFlexmix FLXMRglm
#' @param counts_matrix an \code{\link{alleleCounts}} object with ref and alt slots filled
#' @param sample.id character sample.id to fit model on.
#' @param k vector of mixture components to fit
#' @param coverage_threshold exclude sites with total coverage < coverage_threshold
#' @param niter number of iterations to run
#' @param nrep number of repetitions to run
#' @export 
binommix <- function(counts_matrix, sample.id, k, coverage_threshold = 0L, niter = 1000, nrep = 10) {
    # I/O error handling
    if (!inherits(counts_matrix, "alleleCounts") && length(counts_matrix) != 3) stop("Invalid alleleCounts object")
    stopifnot(is.character(sample.id) && length(sample.id) == 1)
    stopifnot(is.integer(coverage_threshold) && coverage_threshold >= 0L)
    if (any(k < 1L | k > 5L ) ) stop("Number of mixture components must be between 1 and 5")
    
    # data set up
    y <- cbind(counts_matrix$alt[sample.id, ], 
               counts_matrix$ref[sample.id,])
    # filter SNPs that are uninformative for MOI, low coverage
    filter_zeros <- rowSums(y) <= coverage_threshold | is.na(rowSums(y))
    y_obs <- y[!filter_zeros, ]
    
    flexmix::initFlexmix(y_obs ~ 1, 
                         k = k, 
                         model = flexmix::FLXMRglm(y_obs ~ ., family = "binomial"),
                         control = list(iter.max = niter,
                                        minprior = 0),
                         nrep = 10)
}

# mse <- function(counts_matrix, sample.id, fitted_model) {
#    
#}

#' Fit binomial mixture  model in non-overlapping genomic windows
#' 
#' @param counts_matrix an \code{\link{alleleCounts}} object with ref and alt slots filled
#' @param sample.id character sample.id to fit model on.
#' @param window_list list of genomic windows
#' @param bparam options to be passed to \code{\link[BiocParallel]{bpapply}}
#' @return modelWindow object 
#' @importFrom flexmix initFlexmix FLXMRglm
#' @importFrom foreach foreach
#' @importFrom iterators iter
binomMixSlide <- function(counts_matrix, sample.id, window_list, bpparam) {
    # need to re-edit to use BioCParallel
    y <- cbind(counts_matrix$alt[sample.id, ], 
               counts_matrix$ref[sample.id, ])
    
    
    foreach(chr=iter(names(window_list))) %:% 
        foreach(window=unique(window_list[[chr]]$window), .packages = "flexmix") %do% {
            window_y <- y[window_list[[chr]]$variant.id[window_list[[chr]]$window == window], ]
            filter_zeros <- window_y[,1]==0 | window_y[,2] == 0
            y_obs <- y[!filter_zeros, ]
            flexmix::initFlexmix(y_obs ~ 1, 
                                 k = 2, 
                                 model = flexmix::FLXMRglm(y_obs ~ ., family = "binomial"),
                                 control = list(iter.max = 500,
                                                minprior = 0))
        }

}

#' Return estimated model parameters
#' 
#' @importFrom flexmix getModel parameters prior
#' @param model stepFlexmix or flexmix object
#' @param k if stepFlexmix choose model with k components
#' @param criterion if stepFlexmix choose model according to information 
#' criterion
#' @export
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


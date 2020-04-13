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
#' @export 
binommix <- function(counts_matrix, sample.id, k, coverage_threshold = 0L, niter = 1000) {
    # I/O error handling
    if (!inherits(counts_matrix, "alleleCounts") && 
        length(counts_matrix) != 3) stop("Invalid alleleCounts object")
    
    stopifnot(is.character(sample.id) && length(sample.id) == 1)
    if (!(sample.id %in% dimnames(counts_matrix$ref)$sample)) stop("sample.id not found in counts_matrix")
    stopifnot(is.integer(coverage_threshold) && coverage_threshold >= 0L)
    if (any(k < 1L | k > 5L | is.na(k)) ) stop("Number of mixture components must be between 1 and 5")
    
    # data set up
    y <- cbind(counts_matrix$alt[sample.id, ], 
               counts_matrix$ref[sample.id,])
    ds <- counts_matrix$dosage[sample.id, ]
    # filter SNPs that are uninformative for MOI (i.e. non hets), low coverage or missing
    keep_snps <- !is.na(rowSums(y)) & rowSums(y) > coverage_threshold & !is.na(ds) & ds == 1
    y_obs <- y[keep_snps, ]
    baf <- y_obs[,1] / rowSums(y_obs)
    fits <- flexmix::initFlexmix(y_obs ~ 1,
                         k = k, 
                         model = flexmix::FLXMRglm(y_obs ~ 1, family = "binomial"),
                         control = list(iter.max = niter,
                                        minprior = 0),
                         nrep = 5)
    structure(list(fits = fits, baf = baf, sample.id = sample.id), 
              class = "moimix")
}


# plot.moimix <- function(moimix_obj, ...) {
#     
# } 

#' #' Fit binomial mixture  model in non-overlapping genomic windows
#' #' 
#' #' @param counts_matrix an \code{\link{alleleCounts}} object with ref and alt slots filled
#' #' @param sample.id character sample.id to fit model on.
#' #' @param window_list list of genomic windows
#' #' @param bparam options to be passed to \code{\link[BiocParallel]{bpapply}}
#' #' @return modelWindow object 
#' #' @importFrom flexmix initFlexmix FLXMRglm
#' #' @importFrom foreach foreach
#' #' @export
#' binomMixSlide <- function(counts_matrix, sample.id, window_list, bpparam) {
#'     # need to re-edit to use BioCParallel
#'     y <- cbind(counts_matrix$alt[sample.id, ], 
#'                counts_matrix$ref[sample.id, ])
#'     
#'     
#'     foreach(chr=iter(names(window_list))) %:% 
#'         foreach(window=unique(window_list[[chr]]$window), .packages = "flexmix") %do% {
#'             window_y <- y[window_list[[chr]]$variant.id[window_list[[chr]]$window == window], ]
#'             filter_zeros <- window_y[,1]==0 | window_y[,2] == 0
#'             y_obs <- y[!filter_zeros, ]
#'             flexmix::initFlexmix(y_obs ~ 1, 
#'                                  k = 2, 
#'                                  model = flexmix::FLXMRglm(y_obs ~ ., family = "binomial"),
#'                                  control = list(iter.max = 500,
#'                                                 minprior = 0))
#'         }
#' 
#' }




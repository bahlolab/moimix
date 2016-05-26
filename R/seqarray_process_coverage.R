# seqarray_process_coverage.R
# Methods for manipulating read counts  from SeqVarGDS objects

#' Process GATK files
#' @importFrom SeqArray seqGetData
processGATK <- function(gdsfile, ref.allele) {
    # GATK using the AD tag to store read count data
    read_counts <-seqGetData(gdsfile, "annotation/format/AD")$data
    sample.id <- seqGetData(gdsfile, "sample.id")
    variant.id <- seqGetData(gdsfile, "variant.id")
    
    
    # for each variant the AD tag contains the ref and alt counts
    nsnps <- ncol(read_counts)
    ref_index <- seq(1, nsnps, by = 2) # odd index are ref counts
    alt_index <- seq(2, nsnps, by = 2) # even index are alt counts
    
    if (is.null(ref.allele)) {
        counts_matrix <- list(ref = matrix(read_counts[, ref_index], 
                                           nrow = length(sample.id),
                                           ncol = length(ref_index),
                                           dimnames = dimnames(read_counts)),
                              alt = matrix(read_counts[, alt_index],
                                           nrow = length(sample.id),
                                           ncol = length(alt_index),
                                           dimnames = dimnames(read_counts)))
        dimnames(counts_matrix$ref)$sample <-  dimnames(counts_matrix$alt)$sample <- sample.id
        dimnames(counts_matrix$ref)$variant <- dimnames(counts_matrix$alt)$variant <-  variant.id
        return(counts_matrix)
    }
    
    else if (ref.allele == 0L) {
        counts_matrix <- matrix(read_counts[, ref_index], 
                                nrow = length(sample.id),
                                ncol = length(ref_index),
                                dimnames = dimnames(read_counts))
        dimnames(counts_matrix)$sample <- sample.id
        dimnames(counts_matrix)$variant <- variant.id
        return(counts_matrix)
    }
    
    else if (ref.allele == 1L) {
        counts_matrix <- matrix(read_counts[, alt_index], 
                                nrow = length(sample.id),
                                ncol = length(alt_index),
                                dimnames = dimnames(read_counts))
        dimnames(counts_matrix)$sample <- sample.id
        dimnames(counts_matrix)$variant <- variant.id
        return(counts_matrix)
        
    } else {
        stop("Invalid ref.allele argument must be 0 or 1.")
    }
}

#' Process varscan files
#' @importFrom SeqArray seqGetData
processVarscan <- function(gdsfile, ref.allele) {
    # varscan is more sensible
    sample.id <- seqGetData(gdsfile, "sample.id")
    variant.id <- seqGetData(gdsfile, "variant.id")
    if (ref.allele == 0L) {
       counts_matrix <- seqGetData(gdsfile, "annotation/format/RD")$data
    }
    if (ref.allele == 1L) {
       counts_matrix <- seqGetData(gdsfile, "annotation/format/AD")$data
    }
    
    if(is.null(ref.allele)) {
        counts_matrix <- list(ref =seqGetData(gdsfile, "annotation/format/RD")$data,
                        alt =seqGetData(gdsfile, "annotation/format/AD")$data)
        dimnames(counts_matrix$ref)$sample <-  dimnames(counts_matrix$alt)$sample <- sample.id
        dimnames(counts_matrix$ref)$variant.id <- dimnames(counts_matrix$ref)$variant <-  variant.id
        return(counts_matrix)
    }

    dimnames(counts_matrix)$sample <- sample.id
    dimnames(counts_matrix)$variant <- variant.id
    counts_matrix
}

#' Return allele counts from a genofile based on read counts
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @details If \code{ref.allele = NULL}, the function returns a list
#' of matrices corresponding to the read counts of the reference and alternate
#' alleles respectively, while if \code{ref.allele = 0L} or \code{ref.allele = 1L}
#' the function returns a  matrix corresponding the reference or alternate
#' read counts. Currently, the AD tag formats from varscan2re and GATK
#' are supported.
#' @return an alleleCounts matrix object, which is a list containg the 
#' @importFrom SeqArray seqSummary
#' @export 
alleleCounts <- function(gdsfile) {
    # I/O checks
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    vars <-seqSummary(gdsfile, check="none", verbose=FALSE)$format$ID
    
    # GATK will have just AD but no RD
    if("AD" %in% vars && !("RD" %in% vars))  {
        counts_matrix <- processGATK(gdsfile, ref.allele = NULL)
    }
    
    else if("AD" %in% vars && "RD" %in% vars) {
        counts_matrix <- processVarscan(gdsfile, ref.allele = NULL)
    }
    
    else {
        stop("No valid annotation/format tag to compute allele counts matrix")
    }
    counts_matrix$dosage <- getDosage(gdsfile)
    structure(counts_matrix, class="alleleCounts")
}

#' sum number of samples/snps meeting coverage threshold scaled by total number
#' of SNPs/proportions
#' @param x integer vector of read counts
#' @param threshold integer lower bound for coverage
#' @param scale integer to scale results. 
scaledProportion <- function(x, threshold, scale) {
    sum(x > threshold, na.rm = TRUE)/ scale
}

#' SNP coverage rates by sample thresholds
#' 
#' @param coverage_list an \code{\link{alleleCounts}} object
#' @param threshold coverage threshold
#' @description Compute the proportion of samples within each SNP that are covered
#' to at least threshold.
#' @return a vector of n.snps length
#' @export 
coverageBySample <- function(coverage_list, threshold) {
    
    stopifnot(inherits(coverage_list, 'alleleCounts') && length(coverage_list) == 3)
    stopifnot(is.finite(threshold) && is.numeric(threshold))
    total_depth <- coverage_list$ref + coverage_list$alt
    nsamples <- nrow(total_depth)
    apply(total_depth, 2, scaledProportion, threshold = threshold, scale = nsamples)
}

#' Sample coverage rates by SNP thresholds
#'
#' @param coverage_list an \code{\link{alleleCounts}} object
#' @param threshold coverage threshold
#' @description Compute the proportion of SNPs within each sample that are covered
#' to at least threshold.
#' @return a vector of n.samples length
#' @export
coverageBySNP <- function(coverage_list, threshold) {
    stopifnot(inherits(coverage_list, 'alleleCounts') && length(coverage_list) == 3)
    stopifnot(is.finite(threshold) && is.numeric(threshold))
    
    total_depth <- coverage_list$ref + coverage_list$alt    
    nsnps <- ncol(total_depth)
    apply(total_depth, 1, scaledProportion, threshold = threshold, scale = nsnps)
}

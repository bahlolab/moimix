# seqarray_process_coverage.R
# Methods for manipulating read counts  from SeqVarGDS objects

#' Process GATK files
#' @importFrom SeqArray seqGetData
processGATK <- function(gdsfile, ref.allele) {
    # GATK using the AD tag to store read count data
    read_counts <-seqGetData(gdsfile, "annotation/format/AD")$data
    # for each variant the AD tag contains the ref and alt counts
    nsnps <- ncol(read_counts)
    ref_index <- seq(1, nsnps, by = 2) # odd index are ref counts
    alt_index <- seq(2, nsnps, by = 2) # even index are alt counts
    
    if (is.null(ref.allele)) {
        return(list(ref = read_counts[, ref_index],
                    alt = read_counts[, alt_index]))
    }
    
    else if (ref.allele == 0L) {
        return(read_counts[, ref_index])
    }
    
    else if (ref.allele == 1L) {
        return(read_counts[, alt_index])
    }
}

#' Process varscan files
#' @importFrom SeqArray seqGetData
processVarscan <- function(gdsfile, ref.allele) {
    # varscan is more sensible
    if (ref.allele == 0L) {
       seqGetData(gdsfile, "annotation/format/RD")$data
    }
    if (ref.allele == 1L) {
       seqGetData(gdsfile, "annotation/format/AD")$data
    }
    
    if(is.null(ref.allele)) {
        return(list(ref =seqGetData(gdsfile, "annotation/format/RD")$data,
                    alt =seqGetData(gdsfile, "annotation/format/AD")$data))
    }
}

#' Return allele counts from a genofile based on read counts
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param ref.allele \code{NULL} a single numeric value.
#' @details If \code{ref.allele = NULL}, the function returns a list
#' of matrices corresponding to the read counts of the reference and alternate
#' alleles respectively, while if \code{ref.allele = 0L} or \code{ref.allele = 1L}
#' the function returns a  matrix corresponding the reference or alternate
#' read counts. Currently, the variant formats from varscan and GATK
#' are supported.
#' @importFrom SeqArray seqSummary
#' @export 
alleleCounts <- function(gdsfile, ref.allele = NULL) {
    # I/O checks
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(ref.allele) | ref.allele %in% c(0L,1L))
    
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$var.name
    
    # GATK will have just AD but no RD
    if("AD" %in% vars && !("RD" %in% vars))  {
        processGATK(gdsfile, ref.allele)
    }
    
    else if("AD" %in% vars && "RD" %in% vars) {
        processVarscan(gdsfile, ref.allele)
    }
    
    else {
        stop("No valid annotation/format tag to compute allele counts matrix")
    }
    
}

#' sum number of samples/snps meeting coverage threshold scaled by total number
#' of SNPs/proportions
scaledProportion <- function(x, threshold, scale) {
    sum(x > threshold, na.rm = TRUE)/ scale
}

#' SNP coverage rates by sample thresholds
#' 
#' @param coverage_list output from \code{\link{alleleCounts}}
#' @param threshold coverage threshold
#' @description Compute the proportion of samples within each SNP that are covered
#' to at least threshold.
#' @return a vector of n.snps length
#' @export 
coverageBySample <- function(coverage_list, threshold) {
    
    stopifnot(inherits(coverage_list, 'list') && length(coverage_list) == 2)
    stopifnot(is.finite(threshold) && is.numeric(threshold))
    total_depth <- coverage_list$ref + coverage_list$alt
    nsamples <- nrow(total_depth)
    apply(total_depth, 2, scaledProportion, threshold = threshold, scale = nsamples)
}

#' Sample coverage rates by SNP thresholds
#'
#' @param coverage_list output from \code{\link{alleleCounts}}
#' @param threshold coverage threshold
#' @description Compute the proportion of SNPs within each sample that are covered
#' to at least threshold.
#' @return a vector of n.samples length
#' @export
coverageBySNP <- function(coverage_list, threshold) {
    stopifnot(inherits(coverage_list, 'list') && length(coverage_list) == 2)
    stopifnot(is.finite(threshold) && is.numeric(threshold))
    
    total_depth <- coverage_list$ref + coverage_list$alt    
    nsnps <- ncol(total_depth)
    apply(total_depth, 1, scaledProportion, threshold = threshold, scale = nsnps)
}
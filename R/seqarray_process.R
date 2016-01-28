# Title: seqarray_process.R
# Author: Stuart Lee
# Description: Helper functions for processing 
# output from seqarrays GDS format for use in moimix.

processGATK <- function(gdsfile, ref.allele) {
    # GATK using the AD tag to store read count data
    read_counts <- SeqArray::seqGetData(gdsfile, "annotation/format/AD")$data
    # for each variant the AD tag contains the ref and alt counts
    nsnps <- ncol(read_counts)
    ref_index <- seq(1, nsnps*2, by = 2) # odd index are ref counts
    alt_index <- seq(2, nsnps*2, by = 2) # even index are alt counts
    
    if (ref.allele == 0L) {
        return(read_counts[, ref_index])
    }
    
    if (ref.allele == 1L) {
        return(read_counts[, alt_index])
    }
    
    if (is.null(ref.allele)) {
        return(list(ref = read_counts[, ref_index],
                    alt = read_counts[, alt_index]))
    }
    
}

processVarscan <- function(gdsfile, ref.allele) {
    # varscan is more sensible
    if (ref.allele == 0L) {
        SeqArray::seqGetData(gdsfile, "annotation/format/RD")$data
    }
    if (ref.allele == 1L) {
        SeqArray::seqGetData(gdsfile, "annotation/format/AD")$data
    }
    
    if(is.null(ref.allele)) {
        return(list(ref = SeqArray::seqGetData(gdsfile, "annotation/format/RD")$data,
                    alt = SeqArray::seqGetData(gdsfile, "annotation/format/AD")$data))
    }
}

#' Return allele counts from a genofile based on read counts
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param vcf.caller a character string with variant caller name, currently 
#' supports 'gatk_ug', 'gatk_hc', 'varscan'
#' @param ref.allele \code{NULL} a single numeric value.
#' @details If \code{ref.allele = NULL}, the function returns a list
#' of matrices corresponding to the read counts of the reference and alternate
#' alleles respectively, while if \code{ref.allele = 0L} or \code{ref.allele = 1L}
#' the function returns a  matrix corresponding the reference or alternate
#' read counts.
#' @importClassesFrom SeqArray SeqVarGDSClass
#' @export 
alleleCounts <- function(gdsfile, vcf.caller, ref.allele = NULL) {
    # I/O checks
    stopifnot(inherits(genofile, "SeqVarGDSClass"))
    stopifnot(is.null(ref.allele) | ref.allele %in% c(0L,1L))
    
    vars <- seqSummary(gdsfile, "annotation/format", check="none", verbose=FALSE)$ID
    
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

#' Get major allele for isolates
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' 
#' @description Extract major alleles from a GDS file
#' using the Reference allele frequency (RAF) estimated from 
#' read counts directly. In the case of RAF = 0.5 we select the
#' major allele randomly from Bernoulli RV. 
callMajor <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate RAF matrix
    raf_matrix <- 1 - getBAF(gdsfile)
    gt_matrix <- matrix(NA, nrow = nrow(raf_matrix), ncol = ncol(raf_matrix))
    
    gt_matrix
    
}

#' Compute B-allele frequency spectrum
#'
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @import SeqArray
getBAF <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate NRAF matrix, currently on GATK vcf file support
    vars <- seqSummary(gdsfile, "annotation/format", check="none", verbose=FALSE)$ID
    if(!("AD" %in% vars)) {
        stop("Must have annotaion/format/AD tag to compute B-allele frequencies")
    }
    
    seqApply(gdsfile, "annotation/format/AD",
             function(x) x[,2] / rowSums(x),
             margin = "by.variant",
             as.is = "list")
    
}

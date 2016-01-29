# Title: seqarray_process.R
# Author: Stuart Lee
# Description: Helper functions for processing 
# output from seqarrays GDS format for use in moimix.

processGATK <- function(gdsfile, ref.allele) {
    # GATK using the AD tag to store read count data
    read_counts <- SeqArray::seqGetData(gdsfile, "annotation/format/AD")$data
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

#' Get major allele for isolates from BAF matrix
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param get.nucleotides FALSE logical indicating whether to replace allele
#' with nucleotide bases
#' 
#' @description Extract major alleles from a GDS file
#' using the non-reference allele frequency (NRAF) estimated from 
#' read counts directly. In the case of RAF = 0.5 we select the
#' major allele randomly from Bernoulli distribution.
#' @return a binary matrix containing reference and alternate alleles
#' or if get.nucleotides a character matrix
#' @importFrom SeqArray ref alt
#' @export
callMajor <- function(gdsfile, get.nucleotides = FALSE) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate BAF matrix
    baf <- getBAF(gdsfile)
    gt_matrix <- baf
    # if baf is less than 0.5, we'll assign reference allele
    gt_matrix[baf < 0.5] <- 0
    # if baf is greater than 0.5 we'll assign alternate allele
    gt_matrix[baf > 0.5] <- 1
    # if baf is between then we'll assign randomly
    # first need to find bafs == 0.5 and not missing
    index <- baf == 0.5 & !is.na(baf)
    gt_matrix[index] <- rbinom(n = length(gt_matrix[index]), 
                               size = 1, 
                               prob = 0.5)
    
    if (get.nucleotides) {
        # generate character vectors containing nucleotides
        ref.alleles <- as.character(ref(gdsfile))
        alt.alleles <- sapply(1:ncol(gt_matrix), 
                              function(i) as.character(alt(variants)[[i]]))
        # pre-assign matrix to fill
        nt_matrix <- matrix(NA, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
        # array indexes for bases, second column will contain
        # index for character vector assigned above
        ref.indices <- which(gt_matrix == 0, arr.ind = TRUE)
        alt.indices <- which(gt_matrix == 1, arr.ind = TRUE)
        miss.indices <- which(is.na(gt_matrix), arr.ind = TRUE)
        # fill matrix
        nt_matrix[ref.indices] <- ref.alleles[ref.indices[,2]]
        nt_matrix[alt.indices] <- alt.alleles[alt.indices[,2]]
        nt_matrix[miss.indices] <- "N"
        return(nt_matrix)
    }
    
    return(gt_matrix)

}

#' Compute B-allele frequency spectrum
#'
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @importFrom  SeqArray seqSummary seqApply
#' @export
getBAF <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate NRAF matrix, currently on GATK vcf file support
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$var.name
    if(!("AD" %in% vars)) {
        stop("Must have annotaion/format/AD tag to compute B-allele frequencies")
    }
    
    # compute BAF for each sample 
    nrf <- seqApply(gdsfile, "annotation/format/AD",
             function(x) x[,2] / rowSums(x),
             margin = "by.variant",
             as.is = "list")
    # convert list to matrix
    baf <- matrix(unlist(nrf), ncol = length(nrf))
    
    baf
}

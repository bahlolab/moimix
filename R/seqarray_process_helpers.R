# Title: seqarray_process_helpers.R
# Author: Stuart Lee
# Description: Helper functions for processing 
# output from seqarrays GDS format for use in moimix.

#' Helper function for producing data frame of genomic coordinates
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @importFrom SeqArray seqGetData
#' @export
getCoordinates <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    data.frame(chromosome = seqGetData(gdsfile, "chromosome"),
               position = seqGetData(gdsfile, "position"),
               variant.id = seqGetData(gdsfile, "variant.id"),
               stringsAsFactors = FALSE)
}

#' Helper function for producing matrix of allele dosages
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @importFrom SeqArray seqGetData
#' @export
getDosage <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    gt_array <- seqGetData(gdsfile, "genotype")
    sample.id <- seqGetData(gdsfile, "sample.id")
    variant.id <- seqGetData(gdsfile, "variant.id")
    dosage <-  matrix(gt_array[1,,] + gt_array[2,,],
                      nrow = length(sample.id), 
                      ncol = length(variant.id),
                      dimnames = list(sample = sample.id,
                                      variant = variant.id))
    dosage
}

#' Compute per variant depth accross all samples.
#'
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @return an integer vector
perSiteCoverage <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate NRAF matrix, currently on GATK vcf file support
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$ID
    if(!("AD" %in% vars)) {
        stop("gdsfile must have annotation/format/AD tag to compute coverage")
    }
    
    seqApply(gdsfile, "annotation/format/AD", 
             function(x) sum(x, na.rm = TRUE), 
             margin = "by.variant", 
             as.is = "integer")
}

#' Compute per site per sample coverage
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @return a matrix of allele counts
#' @importFrom SeqArray seqSummary seqApply seqGetData
#' @export
perSiteAlleleCoverage <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate NRAF matrix, currently on GATK vcf file support
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$ID
    if(!("AD" %in% vars)) {
        stop("gdsfile must have annotation/format/AD tag to compute ")
    }
    
    matrix(unlist(seqApply(gdsfile, "annotation/format/AD", 
             function(x) colSums(x, na.rm = TRUE), 
             margin = "by.variant", 
             as.is = "list")), ncol = 2, byrow = TRUE)
    
}

#'Compute the proportion of samples covered up any number of bases
#'
#'@param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#'@param threshold integer lower bound for coverage  
#'@return a vector of length equal to number of variants with
#'the proportion of samples covered up to threshold
#'@importFrom SeqArray seqSummary seqApply seqGetData
#'@export
perSiteCoverageBySample <- function(gdsfile, threshold) {
    # i/o checks
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.integer(threshold) && is.finite(threshold))
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$ID
    if(!("AD" %in% vars)) {
        stop("gdsfile must have annotation/format/AD tag to compute ")
    }
    
    # this will produce the depth for each sample at the site
    # i.e. rowSums then sum the number above the threshold
    # then scale by the number of samples
    n.samples <- length(seqGetData(gdsfile, "sample.id"))
    
    seqApply(gdsfile, "annotation/format/AD",
             function(x) sum(rowSums(x) > threshold) / n.samples,
             margin = "by.variant",
             as.is = "double")
    
}

#' Compute proportion of SNPs covered up to a lower bound of bases within 
#' each sample.
#'@param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#'@param threshold integer lower bound for coverage  
#'@return a vector of length equal to number of samples with
#'the proportion of SNPs covered up to threshold
#'@importFrom SeqArray seqSummary seqApply seqGetData
#'@export
perSiteCoverageBySNP <- function(gdsfile, threshold) {
    # i/o checks
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.integer(threshold) && is.finite(threshold))
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$ID
    if(!("AD" %in% vars)) {
        stop("gdsfile must have annotation/format/AD tag to compute ")
    }
    # generate per sample coverage using rowSums(x) for each variant
    # that is greater than threshold 
    # Reduce will over the lists in place then scale  by number of variants.
    n.variants <- length(seqGetData(gdsfile, "variant.id"))
    Reduce("+", 
           seqApply(gdsfile, "annotation/format/AD",
                    function(x) rowSums(x) > threshold,
                    margin = "by.variant",
                    as.is = "list")) / n.variants
    
}

#' Get major allele calls for isolates from BAF matrix
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param get.nucleotides FALSE logical indicating whether to replace allele
#' with nucleotide bases
#' @param use.hets FALSE logical indicating whether to code hets sepeartely
#' 
#' @description Extract major alleles from a GDS file
#' using the non-reference allele frequency (NRAF) estimated from 
#' read counts directly. In the case of RAF = 0.5 we select the
#' major allele randomly from Bernoulli distribution if use.hets is FALSE. Other
#' wise hets are coded as 3. If get.nucleotides is TRUE a character matrix is
#' retruned with the alleles coded, otherwise if use.hets the allele is coded as
#' N. Missing genotype calls are coded as X. 
#' @return a binary matrix  or character matrix containing reference and alternate alleles
#' @importFrom SeqArray ref alt
#' @export
callMajor <- function(gdsfile, get.nucleotides = FALSE, use.hets = FALSE) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate BAF matrix
    baf <- bafMatrix(gdsfile)$baf_matrix
    gt_matrix <- baf
    # if baf is less than 0.5, we'll assign reference allele
    gt_matrix[baf < 0.5] <- 1
    # if baf is greater than 0.5 we'll assign alternate allele
    gt_matrix[baf > 0.5] <- 2
    # if baf is = 0.5 then we'll assign randomly
    # first need to find bafs == 0.5 and not missing
    index <- baf == 0.5 & !is.na(baf)
    if (use.hets) {
        gt_matrix[index] <- 3
    } else {
        gt_matrix[index] <- rbinom(n = length(gt_matrix[index]), 
                                   size = 1, 
                                   prob = 0.5)
    }

    if (!get.nucleotides) {
        return(gt_matrix)
    } else {
        
        # generate character vectors containing nucleotides
        ref.alleles <- as.character(ref(gdsfile))
        alt.alleles <- sapply(1:ncol(gt_matrix), 
                              function(i) as.character(alt(gdsfile)[[i]]))
        # pre-assign matrix to fill
        nt_matrix <- matrix(NA, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
        # array indexes for bases, second column will contain
        # index for character vector assigned above
        ref.indices <- which(gt_matrix == 1, arr.ind = TRUE)
        alt.indices <- which(gt_matrix == 2, arr.ind = TRUE)
        miss.indices <- which(is.na(gt_matrix), arr.ind = TRUE)
        # fill matrix
        nt_matrix[ref.indices] <- ref.alleles[ref.indices[,2]]
        nt_matrix[alt.indices] <- alt.alleles[alt.indices[,2]]
        nt_matrix[miss.indices] <- "X"
        
        if (use.hets) {
            het.indices <- which(gt_matrix == 3, arr.ind = TRUE)
            nt_matrix[het.indices] <- "N"
            return(nt_matrix)
        } else {
            return(nt_matrix)
        }
    }
}

#' Extract PED files from gds
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param moi.estimates a vector of MOI values obtained by \code{\link{binommix}} 
#' for each sample. DEFAULT NULL
#' @param use.hets FALSE include heterozygote genotypes
#' @param outfile prefix of PLINK files for output (default NULL)
#' @details This function writes a plink .ped and .map file for a given
#' gdsfile. If moi.estimates is set then use.hets is redundant. It will set
#' the sex in each the ped file to 2 if MOI > 1 and set heterozygote genotypes
#' to missing for MOI = 1 calls. This is for use in isoRelate.
#' If the use.hets option is true the genotypes are used as is, other
#' wise heterzygotes are set to missing. This function is slow if there
#' are a large number of variants in the GDS file. 
#' @importFrom SeqArray seqGetData seqApply
#' @export
extractPED <- function(gdsfile, moi.estimates = NULL, use.hets = FALSE, outfile = NULL) {
    # construct the ped file requires 6 columns
    # Family ID, Individual ID, Paternal ID, Maternal ID, Sex, Phenotype
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(use.hets) && length(use.hets) == 1)
    if (!(is.null(outfile))) {
        stopifnot(is.character(outfile) && length(outfile == 1))
    }
    
    sample.id <- seqGetData(gdsfile, "sample.id")
    
    if (!is.null(moi.estimates) && length(moi.estimates) != length(sample.id)) {
        stop("Length of moi.estimates must match number of samples")
    }
    
    # extract samples
    # if use.hets we will code hets as missing
    # recode matrix 
    genotypeRecode <- function(gt, hets) {
        if (nrow(gt) != 2) {
            stop("Non-biallelic variant with non-diploid genotype, 
                 please filter before extracting PED file.")
        }
        sites <- t(gt)
        if (hets) {
            return(sites)
        }
        else {
            sites[sites[,1] != sites[,2]] <- NA
            return(sites)
        }
    }
    
    genotypeRecodeMOI <- function(gt, moi) {
        if (nrow(gt) != 2) {
            stop("Non-biallelic variant with non-diploid genotype, 
                 please filter before extracting PED file.")
        }
        
        if(ncol(gt) != length(moi)) {
            stop("Number of samples not equal to number of MOI estimates")
        }
        
        # get dosage per sample
        ds <- colSums(gt)
        # recode het calls as missing for MOI = 1
        recode_hets <- ds == 1 & moi == 1
        
        gt[ , recode_hets] <- NA
        t(gt)
    }
    
    if (!is.null(moi.estimates)) {
        sex <- ifelse(moi.estimates > 1, 2, 1)
        ped_meta <- data.frame(famID = sample.id, indID = sample.id,
                               paternalID = 0, maternalID = 0, 
                               sex = sex, pheno = 2)
        # since we are using moi as sex = 2 we'll code them
        # otherwise for samples where moi = 1 we'll code them as missing
        gt_list <- seqApply(gdsfile, "genotype", 
                            FUN = genotypeRecodeMOI, moi = sex, 
                            margin = "by.variant", as.is = 'list')

        
    } else {
        ped_meta <- data.frame(famID = sample.id, indID = sample.id, 
                               paternalID = 0, maternalID = 0, 
                               sex = 1, pheno = 2)
        gt_list <- seqApply(gdsfile, "genotype", 
                            FUN = genotypeRecode, hets = use.hets, 
                            margin = "by.variant", as.is = 'list')
    }
    
    gt_matrix <- do.call(cbind, gt_list)
    # recode using plink formats
    gt_matrix[gt_matrix == 1] <- 2
    gt_matrix[gt_matrix == 0] <- 1
    gt_matrix[is.na(gt_matrix)] <- 0
    
    ped_data <- cbind(ped_meta, gt_matrix)
    if( !is.null(outfile) ) {
        ped_file <- paste0(outfile, ".ped")
        write.table(ped_data, ped_file, 
                    row.names = FALSE, 
                    col.names = FALSE, 
                    quote = FALSE)
    }

    
    # map file
    chr <- seqGetData(gdsfile, "chromosome")
    pos <- seqGetData(gdsfile, "position")
    snp_id <- paste0(chr, ":", pos)
    # TODO: build proper genetic map for genetic distance setting
    genetic_distance <- pos / 17000
    map_data <- data.frame(chr, snp_id, genetic_distance, pos)
    
    if( !is.null(outfile)) {
        map_file <- paste0(outfile, ".map")
        write.table(map_data, map_file, 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    # return list if the user assigns
    invisible(list(ped = ped_data, map = map_data))
}


#' Get GATK info tags
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' 
#' @description Obtain GATK info tags relevant to best practices.
#' @return a data.frame with the following columns:
#'  chromosome
#'  position
#'  variant.id
#'  qd - (average quality by depth)
#'  mq - RMS mapping quality
#'  fs - Fisher strand bias
#'  mqrank - MQ rank
#'  haplotypescore
#'  readposranksum
#' @importFrom SeqArray seqGetData
#' @references \url{http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set}
#' @export 
getGATKInfo <- function(gdsfile) {
    
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    
    data.frame(chr = seqGetData(gdsfile, "chromosome"),
               pos = seqGetData(gdsfile, "position"),
               variant.id = seqGetData(gdsfile, "variant.id"),
               qd = seqGetData(gdsfile, "annotation/info/QD"),
               mq = seqGetData(gdsfile, "annotation/info/MQ"),
               fs = seqGetData(gdsfile, "annotation/info/FS"),
               mqrank = seqGetData(gdsfile, "annotation/info/MQRankSum"),
               haplotypescore = seqGetData(gdsfile, "annotation/info/HaplotypeScore"),
               readposranksum = seqGetData(gdsfile, "annotation/info/ReadPosRankSum"),
               stringsAsFactors = FALSE)
}

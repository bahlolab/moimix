# Title: seqarray_process_helpers.R
# Author: Stuart Lee
# Description: Helper functions for processing 
# output from seqarrays GDS format for use in moimix.

#' Helper function for producing data frame of genomic coordinates
#' @export
getCoordinates <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    data.frame(chromosome = seqGetData(gdsfile, "chromosome"),
               position = seqGetData(gdsfile, "position"),
               variant.id = seqGetData(gdsfile, "variant.id"),
               stringsAsFactors = FALSE)
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
#' N. 
#' @return a binary matrix  or character matrix containing reference and alternate alleles
#' @importFrom SeqArray ref alt
#' @export
callMajor <- function(gdsfile, get.nucleotides = FALSE, use.hets = FALSE) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate BAF matrix
    baf <- getBAF(gdsfile)
    gt_matrix <- baf
    # if baf is less than 0.5, we'll assign reference allele
    gt_matrix[baf < 0.5] <- 1
    # if baf is greater than 0.5 we'll assign alternate allele
    gt_matrix[baf > 0.5] <- 2
    # if baf is between then we'll assign randomly
    # first need to find bafs == 0.5 and not missing
    index <- baf == 0.5 & !is.na(baf)
    if(use.hets) {
        gt_matrix[index] <- 3
    } else {
        gt_matrix[index] <- rbinom(n = length(gt_matrix[index]), 
                                   size = 1, 
                                   prob = 0.5)
    }
    
    
    if (get.nucleotides) {
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
        }
        
        return(nt_matrix)
    }
    
    return(gt_matrix)

}

#' Extract PED files from gds
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param use.hets FALSE include heterozygote genotypes
#' @param out.file prefix of PLINK files for output
#' @details This function writes a plink .ped and .map file for a given
#' gdsfile. If the use.hets option is true the genotypes are used as is, other
#' wise heterzygotes are set to missing. 
#' @importFrom SeqVarTools getGenotype
#' @export
extractPED <- function(gdsfile, use.hets = FALSE, out.file) {
    # construct the ped file requires 6 columns
    # Family ID, Individual ID, Paternal ID, Maternal ID, Sex, Phenotype
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    sample.id <- seqGetData(gdsfile, "sample.id")
    ped_meta <- data.frame(famID = sample.id, indID = sample.id, 
                           paternalID = 0, maternalID = 0, sex = 1, pheno = 2)
    
    # extract samples
    # if use.hets we will code hets as missing
    # recode matrix 
    genotypeRecode <- function(gt, use.hets) {
        sites <- t(gt)
        if (use.hets) {
            return(sites)
        }
        else {
            sites[sites[,1] != sites[,2]] <- NA
            return(sites)
        }
    }
    # create genotype list
    gt_list <- seqApply(gdsfile, "genotype", 
                        FUN = genotypeRecode, use.hets = use.hets, 
                        margin = "by.variant", as.is = 'list')
    
    gt_matrix <- do.call(cbind, gt_list)
    print(dim(gt_matrix))
    # recode using plink formats
    gt_matrix[gt_matrix == 1] <- 2
    gt_matrix[gt_matrix == 0] <- 1
    gt_matrix[is.na(gt_matrix)] <- 0
    
    ped_file <- paste0(out.file, ".ped")
    ped_data <- cbind(ped_meta, gt_matrix)
    write.table(ped_data, ped_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # map file
    chr <- seqGetData(gdsfile, "chromosome")
    pos <- seqGetData(gdsfile, "position")
    snp_id <- paste0(chr, ":", pos)
    genetic_distance <- pos / 17000
    map_data <- data.frame(chr, snp_id, genetic_distance, pos)
    
    map_file <- paste0(out.file, ".map")
    write.table(map_data, map_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # return list if the user assigns
    invisible(list(ped = ped_data, map = map_data))
}


#' Get GATK info tags
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' 
#' @description Obtain GATK info tags relevant to best practices.
#' @references \url{http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set}
#' @importFrom SeqArray seqGetData
#' @export 
getGATKInfo <- function(gdsfile) {
    
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    
    data.frame(chr = seqGetData(genofile, "chromosome"),
               pos = seqGetData(genofile, "position"),
               variant.id = seqGetData(genofile, "variant.id"),
               qd = seqGetData(genofile, "annotation/info/QD"),
               mq = seqGetData(genofile, "annotation/info/MQ"),
               fs = seqGetData(genofile, "annotation/info/FS"),
               mqrank = seqGetData(genofile, "annotation/info/MQRankSum"),
               haplotypescore = seqGetData(genofile, "annotation/info/HaplotypeScore"),
               readposranksum = seqGetData(genofile, "annotation/info/ReadPosRankSum"),
               stringsAsFactors = FALSE)
}
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
                              function(i) as.character(alt(gdsfile)[[i]]))
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
#' @return a numeric matrix of size l by n where l is the number of samples
#' and n is the number of SNPs. 
#' @export
getBAF <- function(gdsfile, split_by_chrom = FALSE) {
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
    baf <- matrix(unlist(nrf), ncol = length(nrf),
                  dimnames = list(sample = seqGetData(gdsfile, "sample.id"),
                                  variant = seqGetData(gdsfile, "variant.id")))
    
    baf
}

#' Plot B-allele frequencies by sample
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param loic NULL optional list containing chromosome name, start postion, end position, and name
#' @param outdir path to save figures
#' @importFrom scales alpha
#' @export
plotBAF <- function(gdsfile, loci = NULL, outdir) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(dir.exists(outdir))
    if(is.null(loci)) {
        baf <- getBAF(gdsfile)
        # prepare plotting device
        sample.id <- rownames(baf)
        breaks <- tapply(1:ncol(baf), 
                         seqGetData(gdsfile, "chromosome"), median)
        for(sample in sample.id) {
            png(paste0(outdir, "/", sample, "_BAF_all", ".png"))
            plot(baf[sample, ], xaxt ="n", xlab = "", 
                 ylim = c(0,1), ylab = "SNV frequency", 
                 col = scales::alpha("black", 0.5), pch = 16)
            axis(side = 1, at = breaks, labels = names(breaks), 
                 las = 3, cex.axis = 0.6)
            dev.off()
            
        }
    }
    else {
        old.filter <- seqGetFilter(gdsfile)
        seqSetFilterChrom(gdsfile, include = loci$chromosome, is.num = FALSE,
                          from.bp = loci$start, to.bp = loci$end)
        baf <- getBAF(gdsfile)
        sample.id <- rownames(baf)
        
        for(sample in sample.id) {
            png(paste0(outdir, "/", sample, "_BAF_", loci$name, ".png"))
            plot(baf[sample, ], xaxt = "n", xlab = "", ylim = c(0,1), 
                 ylab = "SNV frequency", pch = 16)
            dev.off()
        }
        
        seqSetFilter(gdsfile, variant.sel = old.filter$variant.sel)
        
    }
    
}



#' Compute minor allele frequency 
#' 
#' @details Coverage based estimation of minor allele frequences.
#' The MAF  is computed by computing the minimum of the proportion of 
#' reads covering the reference and alternate allele over all samples.   
#' @note Currently only supports gds files obtained from GATK callers
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @importFrom SeqArray seqApply seqSummary
#' @export
getMAF <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$var.name
    if(!("AD" %in% vars)) {
        stop("Must have annotaion/format/AD tag to compute minor allele frequencies")
    }
    
    maf <- function(x) {

        depth <- sum(x, na.rm = TRUE)
        # frequency of ref and alt reads
        coverage <- colSums(x, na.rm = TRUE)
        min(coverage) / depth
    }
    res <- seqApply(gdsfile, "annotation/format/AD",
                function(x) maf(x),
                margin = "by.variant",
                as.is = "list")
    unlist(res)
}

#' Compute heterozygosity by sample 
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object 
#' @details Follows the method outlined in Manske et al. 2012. Briefly,
#' for each sample, the reference and and alternate allele frequencies are computed
#' as the proportion of reads covering each allele. Then heterozygosity at a SNP
#' is 1 - (raf^2 + aaf^2) 
#' @return a numeric matrix of size l by n where l is the number of samples
#' and n is the number of SNPs.
getHeterozygosityBySample <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate NRAF matrix, currently on GATK vcf file support
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$var.name
    if(!("AD" %in% vars)) {
        stop("Must have annotaion/format/AD tag to compute B-allele frequencies")
    }
    
    # function for computing heterozygosity at variant site
    heterozygosity <- function(x) {
        depth <- rowSums(x)
        
        p <- x[,1] / depth
        q <- x[,2] /depth
        1 - (p^2 + q^2)
        
    }
    het <- seqApply(gdsfile, "annotation/format/AD", heterozygosity, 
             margin = "by.variant", as.is = "list")
    
    het <- matrix(unlist(het), ncol = length(het))
    
    het
}

#' Compute population level heterozygosity
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object 
#' @note Need to update so population ids can be submitted.
getHeterozygosity <- function(gdsfile) {
    # first compute MAF 
    maf <- getMAF(gdsfile)
    1 - (maf^2 + (1-maf)^2)
    
}

#' Compute \eqn{Fws} within-host diversity statistic
#' 
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object 
#' @details Compute the within host diversity statistic according to the
#' method devised in  Manske et.al, 2012. Briefly, within sample heterozygosity
#' and within population heterozygosity are computed and assigned to ten equal sized
#' MAF bins [0.0.05]...[0.45,0.5]. For each bin the mean within sample and population
#' heterozygosity is computed. A regression line of these values through the orgin
#' is computed for each sample. The \eqn{Fws} is then \eqn{1 - \beta}.
#'
#' @references Manske, Magnus, et al. "Analysis of Plasmodium falciparum 
#' diversity in natural infections by deep sequencing." Nature 487.7407 (2012): 
#' 375-379.
#' @note Currently only works on GATK derived gdsfiles. Needs to be updated
#' to define populations.
#' @seealso \code{\link{getHeterozygosity}}, \code{\link{getHeterozygosityBySample}}
#' @export
getFws <- function(gdsfile) {
    # compute heterozygosity at sample level
    sample.het <- getHeterozygosityBySample(gdsfile)
    # compute heterozygosity at population level
    population.het <- getHeterozygosity(gdsfile)
    # create MAF bins, 0..0.05, 0.05..0.1,...,0.45..0.5 
    maf.bins <- findInterval(getMAF(gdsfile), seq(0, 0.5, length.out = 21))
    # find mean population heterozygoisty values
    mu.population.het <- tapply(population.het, bins, mean)
    
    # find mean sample heterozygosity values
    mu.sample.het <- apply(sample.het, 1, 
                           function(x) tapply(x, bins, mean, na.rm = TRUE))
    
    fws <- apply(mu.sample.het, 2, function(x) 1-lm(x ~ mu.population.het -1)$coeff)
                             
    fws
}

#' Estimate variance in BAF spectra along the genome in non-overlapping windows
#' 
#' 
getBAFvar <- function(gdsfile, window.size) {
    # checks
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.numeric(window.size) & length(window.size) == 1)
    stopifnot(is.finite(window.size))
    
    # step 1 -  retrieve BAF matrix
    baf <- getBAF(gdsfile)
    
    # step 2 -  contstruct windows by chromosome
    generateWindows <- function(variant.id, positions, window.size) {
        start <- min(positions)
        end <- max(positions)
        data.frame(variant.id = variant.id, 
                   position = positions,
                   window =findInterval(positions, seq(start, end, by = window.size)))
    }
    coord <- getCoordinates(gdsfile)
    # split by  chromosome
    
    coord_by_chrom <- split(coord, coord$chromosome)
    intervals <- lapply(coord_by_chrom, 
                        function(y) generateWindows(y$variant.id, y$position, window.size))
    
    intervals
    
}



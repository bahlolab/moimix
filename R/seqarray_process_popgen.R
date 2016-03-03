# seqarray_process_popgen.R
# Population genetics methods for SeqVarGDS objects estimated directly from
# read counts

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
#' @importFrom  SeqArray seqSummary seqApply seqGetData
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
    
    het <- matrix(unlist(het), ncol = length(het),
                  dimnames = list(sample = seqGetData(gdsfile, "sample.id"),
                                  variant = seqGetData(gdsfile, "variant.id")))
    
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
    maf.bins <- findInterval(getMAF(gdsfile), seq(0, 0.5, length.out = 11))
    # find mean population heterozygoisty values
    mu.population.het <- tapply(population.het, maf.bins, mean)
    
    # find mean sample heterozygosity values
    mu.sample.het <- apply(sample.het, 1, 
                           function(x) tapply(x, maf.bins, mean, na.rm = TRUE))
    
    fws <- apply(mu.sample.het, 2, function(x) 1-lm(x ~ mu.population.het -1)$coeff)
    
    fws
}
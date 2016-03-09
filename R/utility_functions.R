# Title: utility_functions.R
# Author: Stuart Lee, Katherine Smith
# Description: Utility functions for the haldane package
#--- Datasets
#' Minor Allele Frequenies of ~ 86,000 Plasmodium falciparum SNPs.
#'
#' A dataset containing information about exonic SNPs in Plasmodium falciparum
#' from Manske et al., 2012 \url{http://www.ncbi.nlm.nih.gov/pubmed/22722859}.
#'
#' @name global_pfalciparum_snps
#' @docType data
#' @format A data frame with 86158 rows and 11 variables:
#' \describe{
#'   \item{ID}{Row number}
#'   \item{Chromosome}{Chromosome}
#'   \item{Position}{Physical position of SNP}
#'   \item{snpName}{SNP identifier}
#'   \item{Geme}{Gene name}
#'   \item{MAF_AFR}{Minor allele frequency in Africa}
#'   \item{MAF_SEA}{Minor allele frequency in South East Asia}
#'   \item{MAF_PNG}{Minor allele frequency in Papua New Guinea}
#'   \item{GeneID}{Gene Identifier}
#'   \item{GeneAliases}{Alternate Gene Identifier}
#'   \item{GeneDescription}{Proposed function of gene}
#'
#'
#' }
#' @source \url{http://www.malariagen.net/data}
NULL

#' Plasmodium falciparum SNV data from Gambia.
#' 
#' Variant calls of polymorphic variants from 65 P.falciparum isolates in The Gambia,
#' from Amambua-Ngwa A, Tetteh KK et al., 2012 \url{http://www.ncbi.nlm.nih.gov/pubmed/23133397}.
#' The data was extracted from VCF files and converted to GDS format from Malaria Genomics Consortium Pf3k version 4. See the README file
#' at \url{ftp://ngs.sanger.ac.uk/production/pf3k/release_4/} for details. The GDS object was obtained
#' by filtering the MalariaGen VCF files and keeping only biallelic SNPs that were polymorphic in The Gambia
#' and passed all filters.  
#' 
#' @name gambian_isolates
#' @docType data
#' @format A \code{\link[SeqArray]{SeqVarGDSClass}} object.
#' @source \url{ftp://ngs.sanger.ac.uk/production/pf3k/release_4/}
NULL

#' @useDynLib moimix
#' @importFrom Rcpp sourceCpp
NULL
#> NULL
# Title: coil_functions.R
# Author: Stuart Lee
# Description: Helper functions for extracting barcodes
# for use in COIL - see 
# http://www.broadinstitute.org/infect/malaria/coil/




#' Extract barcode for use in COIL program
#' 
#' 
extractBarcode <- function(gdsfile, variant.id) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # get current filter
    old.filter <- seqGetFilter(gdsfile)
    seqSetFilter(gdsfile, variant.id = variant.id)
    
    majorAlleles <- callMajor(gdsfile)
    
    
    
}
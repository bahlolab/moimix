# seqarray_process_baf_methods.R
# Methods for estimating and plotting BAF spectra

#' Compute B-allele frequency spectrum genome-wide
#'
#' @param gdsfile a \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @importFrom  SeqArray seqSummary seqApply seqGetData
#' @details Constructs a bafMatrix object which is a list consisting of 
#' a matrix of estimated B-allele frequencies and 
#' @return a bafMatrix object 
#' @export
bafMatrix <- function(gdsfile) {
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # estimate NRAF matrix, currently on GATK vcf file support
    vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$var.name
    if(!("AD" %in% vars)) {
        stop("Must have annotation/format/AD tag to compute B-allele frequencies")
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
    # return class of baf
    baf_matrix <- structure(list(baf_matrix = baf, 
                                    coords = getCoordinates(gdsfile)),
                               class = "bafMatrix")
    baf_matrix
}


#' Plot bafMatrix object
#' 
#' @param baf a bafMatrix object
#' @param sample.id character name of sample to plot
#' @details plots the genome-wide signal of MOI within an isolate
#' @importFrom scales alpha
#' @export
plot.bafMatrix <- function(baf, sample.id, ...) {
    if(!(sample.id %in% rownames(baf$baf_matrix))) {
        stop("sample.id not present in bafMatrix object")
    }
    
    breaks <- tapply(1:ncol(baf$baf_matrix),
                     baf$coords$chromosome, median)
    plot(baf$baf_matrix[sample.id, ], xaxt ="n", xlab = "", 
         ylim = c(0,1), ylab = "SNV frequency", 
         col = scales::alpha("black", 0.5), pch = 16, ...)
    axis(side = 1, at = breaks, labels = names(breaks), 
         las = 3, cex.axis = 0.6, ...)
}



generateWindows <- function(variant.id, positions, window_size) {
    start <- min(positions)
    end <- max(positions)
    data.frame(variant.id = variant.id, 
               position = positions,
               window =findInterval(positions, seq(start, end, by = window_size)))
}

averageVar <- function(window, baf_matrix, by.sample) {
    if(length(window$variant.id) > 1 ) {
        sample.var <- apply(baf_matrix[, window$variant.id], 1, 
                            var, na.rm = TRUE)
        if (by.sample) {
            return(sample.var)
        } else {
            mean(sample.var, na.rm = TRUE)
        }
        
    } else {
        NA
    }
}

#' Estimate variance in BAF spectra along the genome in non-overlapping windows
#' 
#' @param baf_matrix a \code{\link[moimix]{bafMatrix}} object
#' @param window_size integer size of window in bp
#' @param by.sample FALSE partition by sample
#' @details This function computes 
#' @return data.frame with chromosome, window id, start, midpoint and end of window
#' and estimates of average variance for window. If by.sample is TRUE, then there
#' will be additional sample.id column with   
#' @export 
getBAFvar <- function(baf_matrix, window_size, by.sample = FALSE) {
    # checks
    stopifnot(inherits(baf_matrix, "bafMatrix"))
    stopifnot(is.numeric(window_size) & length(window_size) == 1)
    stopifnot(is.finite(window_size))
    
    # step 1 -  retrieve BAF matrix
    baf <- baf_matrix$baf_matrix
    
    # step 2 -  contstruct windows by chromosome
    coord <- baf_matrix$coords
    # split by  chromosome
    coord_by_chrom <- split(coord, coord$chromosome)
    intervals <- lapply(coord_by_chrom, 
                        function(y) generateWindows(y$variant.id, y$position, window_size))
    
    # further split list by windows
    intervals_by_window <- lapply(intervals, function(y) split(y, y$window))
    
    # compute the median position for each window for plotting purposes
    median_pos <- lapply(intervals_by_window, 
                         function(chrom) lapply(chrom, 
                                                function(window) data.frame(start = min(window$position),
                                                                            end = max(window$position),
                                                                            mid = median(window$position))))
    median_pos <- do.call(rbind, lapply(median_pos, 
                                        function(x) do.call(rbind, x)))
    
    ids <- matrix(unlist(strsplit(rownames(median_pos), split = "\\.")), 
                  ncol = 2, byrow = TRUE)
    median_pos$chr <- ids[,1]
    median_pos$window <-ids[,2]
    rownames(median_pos) <- NULL
    
    # now apply variance to each window in the list
    baf_var <- lapply(intervals_by_window, 
                      function(chrom) lapply(chrom, 
                                             function(window) averageVar(window, baf, by.sample)))
    if (by.sample) {
        validDF <- function(x) {
            if (length(x) > 1) {
                summary.df <- data.frame(sample.id = names(x), vb = unlist(x))
                rownames(summary.df) <- NULL
                summary.df
            }
        }
        baf_var <- lapply(baf_var, 
                          function(chrom) lapply(chrom, function(x) validDF(x)))
        
        baf_var_df <- do.call(rbind, lapply(baf_var, 
                                            function(y) do.call(rbind, y)))
        ids <- matrix(unlist(strsplit(rownames(baf_var_df), split = "\\.")), 
                      ncol = 3, byrow = TRUE)
        baf_var_df$chr <- ids[,1]
        baf_var_df$window <- ids[,2]
        rownames(baf_var_df) <- NULL
        # merge in 
        return(merge(median_pos, baf_var_df, by = c("chr", "window")))
        
    }
    
    # much simpler if we aren't doing this by sample
    baf_var_df <- do.call(rbind, lapply(baf_var, 
                                        function(x) data.frame(vb = unlist(x))))
    
    baf_var_df$chr <- ids[,1]
    baf_var_df$window <- ids[,2]
    rownames(baf_var_df) <- NULL
    # merge in 
    merge(median_pos, baf_var_df, by = c("chr", "window"))
    
}
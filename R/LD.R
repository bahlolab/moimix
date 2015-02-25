# R Functions for Linkage Disequilibrium
# for Plasmodium
# Stuart Lee
# 01/07/2014

# returns TRUE if any missing values found in vector
any.na <- function(x) {
  any(is.na(x))
}
# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

mergeFiles <- function(fileList) {
  # merge all files in a list into one data frame
  for(file in fileList){
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- read.table(file, header=TRUE, stringsAsFactors = FALSE)
    }
    
    # if the merged dataset does exist, append to it
    if (exists("dataset")){
      temp_dataset <-read.table(file, header=TRUE, stringsAsFactors = FALSE)
      dataset<-rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
  }
  return(dataset)
}
# Assuming the dataset is a matrix consisting of SNP ids as columns
# and rows are represented by samples

getChrSNPs <- function(data, pattern) {
  # return a vector of SNP ids corresponding to a chromosome
  ids <- colnames(data)[2:ncol(data)]
  
  #pattern <- paste0("_", sprintf("%02d", chrNo), "_")
  
  ids[grepl(pattern, ids)]
}

getChrIndex <- function(data, pattern) {
  colIndex <- grep(pattern, colnames(data))
  colIndex
}

countSNPsChr <- function(data, pattern) {
  # count the number of SNPs in a chromosome
  sum(grepl(pattern, colnames(data)[2:ncol(data)]))
}

missSNP <- function(data, snp) {
  # returns a boolean if any values in the sample 
  # are missing for a SNP
  any(is.na(data[,snp]))
}

nmissSNP <- function(data, snp) {
  # count the number of missing values on a SNP id
  sum(is.na(data[,snp]))
}

bpSNP <- function(snpID) {
  # the base position is at the end of the string
  match <- regexpr("[[:digit:]]{3,}", snpID)
  match.length <- attr(match, "match.length")
  substr(snpID, match, match +match.length)
}

chrSNP <- function(snpID) {
  # chromosome position is 
  match <- regexpr("_[[:digit:]]{2}", snpID)
  match.length <- attr(match, "match.length")
  substr(snpID, match + 1, match + match.length -1)
}

maf <- function(data, snp) {
  # find minor allele frequency of SNP
  n <- length(data[, snp]) - sum(is.na(data[,snp]))
  a1 <- sum(data[, snp] == 1, na.rm = TRUE)/n
  a2 <- sum(data[, snp] == 2, na.rm = TRUE)/n
  return(min(a1, a2))
}

positionSNPs <- function(data, chrNo) {
  # create a data.frame with all snp IDs, base position, chromosome number
  lookupTable <- data.frame(snpIDs = getChrSNPs(data, paste0("_", sprintf("%02d", chrNo), "_")), 
                            stringsAsFactors = FALSE)
  lookupTable$bp <- as.numeric(bpSNP(lookupTable$snpIDs))
  return(lookupTable[order(lookupTable$bp),])
}

chisqPval <- function(data, snp1, snp2) {
  #levels <- c(1,2)
  #x <- factor(data[,snp1], levels = levels)
  #y <- factor(data[,snp2], levels = levels)
  tab <- table(data[, snp1], data[,snp2])
  suppressWarnings(chisq.test(tab)$p.value)
}

chisqPvalAdjusted <- function(data, snp1, snp2, mc = "holm", ncomp) {
  # find the pairwise p.value between two SNPs
  # data: data frame with SNPs
  # snp1, snp2: string identifer or colno for SNP
  # adjust for multiple comparison, default method is holm
  suppressWarnings(p.adjust(chisq.test(data[,snp1], data[,snp2])$p.value, 
                            method = mc, n = ncomp))
}


allLD <- function(data, chrNo, FUN = LD, ...) {
  # find all the pairwise LD values accross a chromosome
  # output a matrix of SNP pairs
  
  colIndex <- getChrIndex(data, paste0("_", sprintf("%02d", chrNo), "_"))
  #chrData <- data[ , getChrSNPs(data, paste0("_", sprintf("%02d", chrNo), "_"))]
  ld <- function(snp1, snp2){
    FUN(data, snp1, snp2, ...)
  }
  snp <- colnames(data)[colIndex]
  results <- outer(colIndex, colIndex, FUN = Vectorize(ld))
  rownames(results) <- snp
  colnames(results) <- snp
  #results[upper.tri(results)] <- 0
  return(results)
}

heatmapLD <- function(mat, ...) {
  # Create a heatmap from a matrix of LD values
  require(LDheatmap)
  distances <- as.numeric(bpSNP(colnames(mat)))
  LDheatmap(mat, genetic.distances = distances, 
            SNP.name = rownames(mat), ...)
}

LD <- function(data, snp1, snp2, type = "r2") {
  # find the linkage disequilbirum between two SNPs
  # return
  if(!(type %in% c("r2", "D`", "chisq"))) {
    stop(paste(type, "is not a valid LD measure."))
  }
  # tabulate  and transform into frequencies
  x <- factor(data[, snp1], levels = c(1,2))
  y <- factor(data[, snp2], levels = c(1,2))
  cTab <- table(x, y)
  #print(c(snp1, snp2, dim(cTab)))
  
  # get p-value
  if(type == "chisq") {
    return(suppressWarnings(chisq.test(cTab)$p.value))
  }
  
  total <- sum(cTab)
  cTab <- cTab / total
  
  # find the LD coeffecient
  D <- cTab[1,1]*cTab[2,2] - cTab[1,2]*cTab[2,1]
  #print(D)
  
  # find allele frequencies
  pA <- sum(cTab[1,])
  qa <- 1 - pA
  
  pB <- sum(cTab[,1]) 
  qb <- 1 - pB
  #print(cbind(pA,qa, pB, qb))
  if(type == "D`") {
    #print(D)
    if(D > 0) {
      scale <- min(pA*qb, qa*pB)
      #print(scale)
    }
    else if(D < 0){
      scale <- max(-pA*pB, -qa*qb)
    }
    
    else {
      return(D)
    }
    
    return(D / scale)
  }
  
  if(type == "r2") {
    return(D^2 / (pA*qa*pB*qb))
  }
  
}

testLD <- function(data, test.data, program = "PLINK", type = "r2") {
  # test results of plink/haploview output for ld measure
  # against plink output
  # data: data frame of SNPs
  # test.data: results from plink, or haploview LD output file
  if(program == "PLINK") {
    # PLINK only outputs r2 genomewide
    testSNPs <- test.data[, c("SNP_A", "SNP_B", "R2")]
    colnames(testSNPs)[3] <- type
    # plink outputs results to 6 decimal
    tolerance <- 1e-6
    #print(head(testSNPs))
  }
  
  if(program == "haploview") {
    # haploview can output either r2 or D`
    # L1 and L2 columns are the SNP names
    testSNPs <- test.data[, c("L1", "L2", "D.", "r.2")]
    colnames(testSNPs) <- c("SNP_A", "SNP_B", "D`", "r2")
    # haploview output is rounded to 3 decimal places
    tolerance <- 1e-3
    
    #print(head(testSNPs))
  }
  #for each pair in testSNPs find LD between SNPs
  print(head(testSNPs))
  testSNPs$results.LD <- NULL
  
  for(i in 1:nrow(testSNPs)){
    testSNPs$results.LD[i] <- LD(data, testSNPs[i,"SNP_A"], testSNPs[i,"SNP_B"], type = type)
  }
  
  
  print(paste("A selection of ", program, " results vs. my LD function."))
  print(testSNPs[sample(1:nrow(testSNPs), 30),])
  all.equal(testSNPs[, type], testSNPs$results.LD, tolerance = tolerance)
}

findPairs <- function(obj) {
  # find pairs of object elements
  res = list()
  for(i in 1:length(obj)){
    for(j in i:length(obj)){
      if(obj[i] != obj[j]) {
        res[[length(res)+1]] <- c(obj[i], obj[j])
      }
    }
  }
  return(res)
}

multiLD <- function(data, snp1, snp2) {
  c(r2 = LD(data, snp1, snp2), dprime = LD(data, snp1, snp2, type="D`"), 
    pvalue = chisqPval(data, snp1, snp2))
}

windowedLD <- function(data, chrNo, chrSize, window) {
  lookupTable <- positionSNPs(data, chrNo)
  regions <- seq(1, chrSize, window)
  for(i in 1:(length(regions) -1)) {
    # find SNPs that lie in the window region
    ranges <- regions[i] <= lookupTable$bp & lookupTable$bp <= regions[i+1]
    snpWindow <- lookupTable$snpIDs[ranges]
    mp <- (regions[i] + regions[i+1])/2
    
  }
  
}

findLDChr <- function(data, chrNo, windowSize = 1000, step = 500) {
  # naive method for finding regions of LD by shuffling along a window
  # and averaging over that region
  # At the moment, the windowing measure doesn't take into account
  # the variance produced by small numbers of SNPs in a window
  # the variance produced by the window size and step size
  # Inputs:
  #   data: set of SNP markers
  #   chrNo: chromosome of interest
  #   windowSize: window size (physical distance)
  #   step: step size (physical distance)
  #   FUN: LD function to calculate
  # Outputs:
  #   results: an n by 4  numeric matrix,
  #   the first column is the midpoint of the window
  #   the following columsn are the LD measures
  lookupTable <- positionSNPs(data, chrNo)
  #print(lookupTable)
  
  # initialise windows
  startPos <- min(lookupTable$bp)
  stopPos <- max(lookupTable$bp) - windowSize  
  windows <- seq(from = startPos, to = stopPos, by = step)
  
  # initialise results matrix
  n <- length(windows)
  results <- matrix(nrow = n, ncol = 4)
  j <- 1
  
  # LD measure
  S <- function(snp1, snp2){
    multiLD(data, snp1, snp2)
  }
  
  # loop over windows
  for(i in windows){
    #startPosition <- lookupTable$bp[i]
    # find midpoint of window
    mp <- i + windowSize/2
    # find index of SNPs that are within the window size
    ranges <- (i <= lookupTable$bp) & (lookupTable$bp <= (i + windowSize))
    
    # if no other SNPs are in the window, continue
    snpWindow <- lookupTable$snpIDs[ranges]
    if(length(snpWindow) < 2) {
      next
    }
    #print(c(i, mp, i+windowSize))
    
    # number of SNPs in the window
    nWindow <- length(snpWindow)
    # all pairwise combinations of SNPs in the window
    snpPairs <- findPairs(snpWindow)
    #print(snpPairs)
    # all pairwise LD of SNPs in the window
    ldPairs <- sapply(snpPairs, FUN = function(x) S(x[1], x[2]))
    # average LD over the window
    Swindow <- apply(ldPairs, 1, sum) / (nWindow * (nWindow -1) / 2)

    # append to results matrix
    results[j,] <- c(mp, Swindow)
    # move to next row of results
    j <- j + 1
  }
  
  # if there were less than 2 SNPs in a window
  # the results matrix was not filled.
  results <- na.omit(results)
  colnames(results) <- c("bp", "r2", "D`", "p.value")
  return(na.omit(results))
}


plotLD <- function(ldmat, chrNo, ...) {
  # Function for producing plots for each LD measure
  # Inputs:
  # ldmat: a matrix of LD measures on a chromsome
  # the chromsome where the LD results were obtained
  
  # sub function for controlling output
  plotWindow <- function(x,y,...) {
    x <- x / 1000
    plot(x, y, xlab="Position (kb)", ylim = c(0,1), cex = 0.6, pch = 4,
         xlim = c(min(x), max(x)), cex.axis = 0.8, cex.lab = 0.8, ...)
    rug(x, ticksize = 0.03, side = 1)
  }
  # r^2 plot
  plotWindow(ldmat[,1], ldmat[,2], 
             cex.main = 0.8, ylab = expression(r^2), 
             main = expression("LD measure: "~ r^2))
  
  # D' plot
  plotWindow(ldmat[,1], ldmat[,3], cex.main = 0.8, ylab = "D`",
             main = expression("LD measure: D`"))
  
  # Raw p-value plot
  plotWindow(ldmat[,1], ldmat[,4], cex.main = 0.8, 
             ylab = expression(Chi[1]^2), 
             main = expression("LD measure: " ~ Chi[1]^2 ~ "raw p-value"))
  #dev.off()
}

writePlink <- function(data, fname, chrNo, haploView = TRUE) {
  # utility function to write .map and .ped file for plink input
  # also write .info file for haploview
  # data: data frame containing SNPs
  # fname : output file name
  # chrNo : chromosome number
  # haploView : boolean (if true write a .info file)
  
  # generate a list of SNPs on the chromosome
  pattern <- paste0("_", sprintf("%02d", chrNo), "_")
  snpIDs <- getChrSNPs(data, pattern)
  
  # find base pair positions
  bpIDs <- bpSNP(snpIDs)
  
  # construct the ped file requires 6 columns
  # Family ID, Individual ID, Paternal ID, Maternal ID, Sex
  # Phenotype
  ped <- data.frame(famID = 1 : nrow(data), indID = 1:nrow(data), 
                    paternalID = 0, maternalID = 0, sex = 1,
                    pheno = 2)
  
  # the ped file also requires genotypes, so for our parasites
  # we need two copies of SNP data and to reset missingness to 0
  snp <- data[, sort(rep(snpIDs,2))]
  snp[is.na(snp)] <- 0
  
  ped <- cbind(ped, snp)
  
  ped.file <- file(paste0(fname, ".ped"), open ="w")

  write.table(ped, file = ped.file, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)  
  
  close(ped.file)
  
  # construct the map file which requires 4 columns as input
  # chromosome (Y because parasites are haploid)
  # snpID
  # distance which we set to 0
  # base pair position, which is given 
  map <- data.frame(chr = "Y", snpIDs = snpIDs, distance = 0, bp = bpSNP(snpIDs))
  map <- map[order(snpIDs),]

  map.file <- file(paste0(fname, ".map"), open = "w")
  write.table(map, map.file, row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
  close(map.file)
  
  if (haploView) {
    
    info.file <- file(paste0(fname, ".info"), open = "w")
    
    write.table(map[, c("snpIDs", "bp")], file = info.file, row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    close(info.file)
  }
}

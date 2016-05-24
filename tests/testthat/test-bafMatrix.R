# test baf-methods
context("bafMatrix object tests")

library(SeqArray)

gds_file <- system.file("extdata", "gambian_isolates.gds", package = "moimix")
isolates <- seqOpen(gds_file, allow.duplicate = TRUE)

# filter to one sample and randomly sample 10000 variants
set.seed(20160523)
test_sample <- "PA0021-C"
n.variants <- length(seqGetData(isolates, "variant.id"))
seqSetFilter(isolates, sample.id = "PA0021-C", 
             variant.id = sample.int(n.variants, size = 10000))



test_that("I/O checks for bafMatrix objects", {
    expect_error(bafMatrix(matrix(rpois(10), nrow = 2)))
    expect_error(bafMatrix("isolates"))
    
})

test_baf <- bafMatrix(isolates)
# for a one sample bafMatrix the baf_matrix should equal baf_site
test_that("bafMatrix sanity checks", {
    expect_equal(as.numeric(test_baf$baf_matrix), test_baf$baf_site)
})

seqClose(isolates)
          
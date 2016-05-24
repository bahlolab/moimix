# test-binomix.R
context("Fitting moimix objects\n")

library(SeqArray)

gds_file <- system.file("extdata", "gambian_isolates.gds", package = "moimix")
isolates <- seqOpen(gds_file, allow.duplicate = TRUE)

# filter to one sample and randomly sample 10000 variants
set.seed(20160524)
n.variants <- length(seqGetData(isolates, "variant.id"))
seqSetFilter(isolates, sample.id = "PA0022-C", 
             variant.id = sample.int(n.variants, size = 10000))

test_sample <- "PA0022-C"
counts_matrix <- alleleCounts(isolates)


test_that("I/O checks for binommix",
          {
              expect_error(binommix(counts_matrix, test_sample, k = -1),
                           "Number of mixture components must be between 1 and 5")
              expect_error(binommix(counts_matrix, test_sample, k = c(NA, 5)),
                           "Number of mixture components must be between 1 and 5")
              expect_error(binommix(list(counts_matrix$ref, counts_matrix$alt), test_sample, k = 3),
              "Invalid alleleCounts object")
              expect_error(binommix(counts_matrix, test_sample, k = 2, coverage_threshold = -20))
              expect_error(binommix(counts_matrix, "abcd", k = 3), "sample.id not found in counts_matrix")
          }
)

seqClose(isolates)
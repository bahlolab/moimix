# Playing with seqArray functions
# work in progress
#this will get genotypes
# seqSlidingWindow(varscan.vcf, "genotype", win.size = 3, shift = 1,
#                  function(x, index) {
#                                      print(x)},
#                  as.is = "none")
#
# seqSlidingWindow(varscan.vcf, "genotype", win.size=4,
#                  FUN = function(index, x) {
#                    z <- unlist(lapply(x, function(z) mean(z, na.rm=TRUE)))
#                    cat("Window ", index, ", starting from Variant ", index,
#                        "\n    ", format(round(z,3), nsmall=3, width=8), "\n", sep="")
#                  },
#                  as.is="none",
#                  var.index="relative")


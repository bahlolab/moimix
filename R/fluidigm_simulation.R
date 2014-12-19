# Title: fluidgim_simulation.R
# Description: Methods and functions for simulating fluidgm assay data
# Details: Two models for simulating data:
#   1. Simulating on a bivariate circular distribution
#   2. Jittering along an arc
# Other ideas:
# Transform space into polar or spherical coordinates?
# Date: 19/14/12


simFluidigm <- function(n.samples, p.ref, ref.allele.info, alt.allele.info) {
  if ( length(ref.allele.info) != 2 | length(alt.allele.info) != 2 ) {
    stop("Reference or alternate allele signal information must be numeric
         vector of lenght 2")
  }
  if ( !is.numeric(p.ref) | (p.ref <= 0 | p.ref >= 1) ) {
    stop("Reference allele probability must be between 0 and 1")
  }

  p.alt <- 1 - p.ref
  p <- sort(c(p.alt, p.ref))
  print(p)
  # generate arc length
  t <- seq(0, 1, length.out = n.samples)
  print(t)
  #runif(n.samples, min = 0, max = p.ref)

  # noise in ref allele
  beta <- rnorm(n.samples,
                mean = ref.allele.info[1],
                sd = ref.allele.info[2])

  # noise in alt allele
  alpha <- rnorm(n.samples,
                 mean = alt.allele.info[1],
                 sd = alt.allele.info[2])
  a <- 0.6
  b <- 0.4
  arc <- cbind(p.alt*sin(alpha*t), p.ref*cos(beta*t))

  arc
}

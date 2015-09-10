# simulation_utility.R
# Description: utility functions to run simulations

#--- Studies
# study 1, no error, fixed coverage of 50, fixed. no.snps, 
# fixed no.samples moi ranging from 2:5,
# random genotype proportions, random clonal proportions
study1 <- function(n.samples, n.snps, ...) {
    out <- list()
    coverage <- rep(50, n.samples)
    out <- foreach(moi = 2:5, .packages = "moimix") %dopar% {
        simulate_moi(n.samples, n.snps, moi, coverage, error = 0, ...)
    }
    
    return(out)
    
}

# study2, same as study 1 with error rate of 1%
# want to make sure you use same settings as study 1,
# so pass study 1 as an input
study2 <- function(n.samples, n.snps, study1.out, ...) {
    out <- list()
    coverage <- rep(50, n.samples)
    out <- foreach(moi = 2:5, .packages = "moimix") %dopar% {
        simulate_moi(n.samples, n.snps, moi, coverage, error = 0.01, 
                     pi.true = study1.out[[moi-1]][["pi.true"]],
                     mu.true = study1.out[[moi-1]][["mu.true"]],
                     aaf = study1.out[[moi-1]][["aaf"]], ...)
    }
}

# study3, same as study 2 but with varying coverage
study3 <- function(n.samples, n.snps, coverage, study2.out, ...) {
    out <- list()
    out <- foreach(moi = 2:5, .packages = "moimix") %dopar% {
        simulate_moi(n.samples, n.snps, moi, coverage = coverage, error = 0.01,
                     pi.true = study2.out[[moi-1]][["pi.true"]],
                     mu.true = study2.out[[moi-1]][["mu.true"]],
                     aaf = study2.out[[moi-1]][["aaf"]], ...)
    }
}

# study4 is slightly different, alter SNP number
# pi is now fixed, mu can varies
study4 <- function(n.snps = c(1000, 10000, 25000, 50000, 100000), ...) {
    pi.true <- matrix(c(0.05, 0.05, 0.9,
                        0.1, 0.2, 0.7,
                        1/3, 1/3, 1/3,
                        0.25, 0.6, 0.15,
                        0.75, 0.2, 0.05,
                        0.66, 0.33, 0.01,
                        0.2, 0.2, 0.6,
                        0.4, 0.4, 0.2,
                        0.24, 0.36, 0.4,
                        0.8, 0.1, 0.1,
                        0.02, 0.03, 0.95,
                        0.5, 0.25, 0.25,
                        0.85, 0.14, 0.01,
                        0.5, 0.4, 0.1,
                        0.55, 0.44, 0.01,
                        0.005, 0.305, 0.69,
                        0.2,0.5,0.3,
                        0.35, 0.25, 0.4,
                        0.99, 0.005, 0.005,
                        0.4, 0.15, 0.45), nrow = 3)
    
    mu.true <- MCMCpack::rdirichlet(20, alpha = rep(1, 3)))
    
    out <- list()
    coverage = rep(50, 20)
    out <- foreach(s = n.snps, .packages = "moimix") %dopar% {
        simulate_moi(n.samples = 20, n.snps = s, moi = 3, coverage = coverage,
                     error = 0.01, pi.true = pi.true, mu.true = mu.true, ...)
    }
    return(out)   
}



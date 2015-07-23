# test-emsteps.R
context("Estep and Mstep testing for binomial mixture models")

test_that(desc = "Mstep improves with more n", {
    # generate sample for mixture distribution
    u = sampleMM(100, 1, 2, mu = c(0.5, 0.6), pi = c(0.9, 0.1))
    # compute E-step to generate sufficient statistics
    # estimate components directly from data
    mixture.weights <- tabulate(u$states) / 100
    mixture.comp <- c(mean(u$obs[u$states==1]), mean(u$obs[u$states==2]))
    ss <- updateClassProb(x = u$obs, N=1, k=2, mixture.comp, mixture.weights)
    theta.small <- mStep(u$obs, 1, ss$tau)
    
    u = sampleMM(1000, 1, 2, mu = c(0.5, 0.6), pi = c(0.9, 0.1))
    mixture.weights <- tabulate(u$states) / 1000
    mixture.comp <- c(mean(u$obs[u$states==1]), mean(u$obs[u$states==2]))
    ss <- updateClassProb(x = u$obs, N=1, k=2, mixture.comp, mixture.weights)
    theta.large <- mStep(u$obs, 1, ss$tau)
    
    dist1.mu <- sum((c(0.5, 0.6) - theta.small$mixture.comp)^2)
    dist2.mu <- sum((c(0.5, 0.6) - theta.large$mixture.comp)^2)
    dist1.pi <- sum((c(0.9, 0.1) - theta.small$mixture.weights)^2)
    dist2.pi <- sum((c(0.9, 0.1) - theta.large$mixture.weights)^2)
    expect_less_than(dist2.mu, dist1.mu)
    expect_less_than(dist2.pi, dist1.pi)
    
})

test_that(desc = "Log-likelihood increases ever iter", {
    
})

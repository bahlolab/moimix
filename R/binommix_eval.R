# binommix_eval.R
# Description: Evaluate the mixture model
# accross a whole bunch of conditions
# Author: Stuart Lee
moi <- 2:5
# seed for MC simulation
set.seed(666351)
B <- 1000 # no. monte-carlo runs
n.snps <- 10000 # no. snps

sim1 <- list()

for(i in moi) {
  sim1[[paste0("moi", i)]] <- simulateMOI(n.samples = B,
                                          n.snps = n.snps,
                                          moi = i,
                                          coverage = rep(50, B),
                                          error = 0,
                                          dirichlet.param = c(1,2),
                                          maf.dist = rbeta,
                                          shape1 = 2,
                                          shape2 = 6)
}


fitEM <- function(sim.list, moi, B) {
  x <- sim.list[[paste0("moi", moi)]][["obs.alt.counts"]]
  true.props <- lapply(sim.list[[paste0("moi", moi)]][["clone.props"]], sort)
  # convert to a matrix
  true.props <- matrix(unlist(true.props), nrow = B, ncol = moi, byrow = TRUE)
  # fit model by ro


  em.fit <- list()

  for (i in 1:B) {
    em.fit[[i]] <- sort(binommixEM(x[i, ], N = 50, k = moi)$pi)
  }

  estimated.props <- matrix(unlist(em.fit),
                            ncol = moi,
                            nrow = B,
                            byrow = TRUE)

  return(list(est.props = estimated.props ,
              true.props = true.props))

}

em.fit <- list()

for (i in moi) {
  em.fit[[paste0("moi", i)]] <- fitEM(sim1, i, B)
}

#-- simulation 2, error in base calling
sim2 <- list()

for (i in moi) {
  sim2[[paste0("moi", i)]] <- simulateMOI(n.samples = B,
                                         n.snps = n.snps,
                                         moi = i,
                                         coverage = rep(50, B),
                                         error = 0.01,
                                         dirichlet.param = c(1,2),
                                         maf.dist = rbeta,
                                         shape1 = 2,
                                         shape2 = 6)
}

em.fit2 <- list()

for (i in moi) {
  em.fit2[[paste0("moi", i)]] <- fitEM(sim2, i, B)
}

#-- evaluation

mseMC <- function(true.mat, est.mat) {
  rowSums((true.mat - est.mat)^2)
}

sim1.results <- lapply(em.fit, function(x) mseMC(x$true.props, x$est.props))
sim2.results <- lapply(em.fit2, function(x) mseMC(x$true.props, x$est.props))

library(reshape2)
sim1.results.melt <- melt(sim1.results)
sim1.results.melt$L1 <- as.numeric(gsub("moi ", " ", sim1.results.melt$L1))
sim1.results.melt$study <- 'no error'
sim2.results.melt <- melt(sim2.results)
sim2.results.melt$L1 <- as.numeric(gsub("moi", " ", sim1.results.melt$L1))
sim2.results.melt$study <- 'error'
# plot results
library(ggplot2)
library(dplyr)
results.all <- rbind(sim1.results.melt, sim2.results.melt)

# produce box plot
p1 <- ggplot(results.all, aes(y = value, x = factor(L1), colour = factor(L1))) +
  geom_boxplot() +
  facet_grid(study ~ .) +
  xlab("") +
  ylab("Mean square error") + scale_colour_discrete(guide = FALSE)

results.all %>% group_by(study, L1) %>% summarise(mse = mean(value), var = var(value) / B,
                                                  lower = mse - 1.96 * sqrt(var),
                                                  upper = mse + 1.96 * sqrt(var)) %>%
  select(L1, mse, lower, upper)

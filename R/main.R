# Bayesian Methods I Project

# Mixture model of two component univariate mixture of normals

# libraries
library(Rcpp)

# simulated data
y <- 0.3*rnorm(100, 0, 1) + 0.7*rnorm(100, 3, 2)

# priors
# thetas ~ N(mu, tau^2)
# precisions ~ gamma(v/2, v*sigma^2/2)
# p ~ beta(a, b)

# theta priors
mu.0 <- 0
tau.20 <- 10

# s priors
sigma.20 <- 1
v.0 <- 5

# p priors
a <- 1
b <- 3

# load C++ MCMC sampler
sourceCpp("src/sampler.cpp")

# number of iterations
S <- 1e5

# explore joint posterior
PHI <- sampler(y, mu.0, tau.20, sigma.20, v.0, a, b, S)

# drop burnins
PHI <- PHI[1001:nrow(PHI), ]

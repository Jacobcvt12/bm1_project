# Bayesian Methods I Project

# Mixture model of two component univariate mixture of normals

# libraries
library(Rcpp)
library(ggplot2)

# read in data
y <- scan("data/simulated.dat")

# priors
# thetas ~ N(mu, tau^2)
# precisions ~ gamma(v/2, v*sigma^2/2)
# p ~ beta(a, b)

# load C++ MCMC sampler
sourceCpp("src/sampler.cpp")

# number of iterations
S <- 1e6

# explore joint posterior
PHI <- sampler(y, mu.1, mu.2, tau.1.2, tau.2.2, a, b, S)

# MCMC diagnostics

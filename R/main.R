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

# theta priors
mu.0 <- 0
tau.20 <- 10

# s priors
sigma.20 <- 1
v.0 <- 5

# p priors
a <- 1
b <- 1

# load C++ MCMC sampler
sourceCpp("src/sampler.cpp")

# number of iterations
S <- 10

# explore joint posterior
set.seed(42)
PHI <- sampler(y, mu.0, tau.20, sigma.20, v.0, a, b, S)

# MCMC diagnostics

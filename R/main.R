# Bayesian Methods I Project

# Mixture model of two component univariate mixture of normals

# read in data
#y <- scan("data/simulated.dat")

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
#sourceCpp("src/sampler.cpp")

# number of iterations
S <- 1e5
B <- 1e4

# explore joint posterior
set.seed(42)
y <- 0.3*rnorm(100)+0.7*rnorm(100, 2, 3)
PHI <- sampler(y, mu.0, tau.20, sigma.20, v.0, a, b, S, B)

# MCMC diagnostics

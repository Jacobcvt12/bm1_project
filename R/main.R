# Bayesian Methods I Project

# Mixture model of two component univariate mixture of normals

# read in data
#y <- scan("data/simulated.dat")
simulateZ <- function(N, p){
    P <- rdirichlet(N, p)
    cumP <- t(apply(P, 1, cumsum))
    u <- runif(N)
    zz <- rep(NA, N)
    zz[u < cumP[, 1]] <- 1
    k <- 2
    while(k <= ncol(P)){
        zz[u < cumP[, k] & u >= cumP[, k-1]] <- k
        k <- k+1
    }
    zz
}
library(gtools)
p <- c(0.3, 0.7)
theta <- c(0, 2)
sigma <- c(1, 3)
z <- simulateZ(100, p)
y <- rnorm(100, theta[z], sigma[z])

# priors
# thetas ~ N(mu, tau^2)
# precisions ~ gamma(v/2, v*sigma^2/2)
# p ~ beta(a, b)

# theta priors
mu.0 <- 0
tau.20 <- 5

# s priors
sigma.20 <- 1
v.0 <- 5

# p priors
a <- 2
b <- 1

# load C++ MCMC sampler
#sourceCpp("src/sampler.cpp")

# number of iterations
S <- 1e5
B <- 1e3

# explore joint posterior
set.seed(42)
#library(coda)
#PHI <- sampler(y, mu.0, tau.20, sigma.20, v.0, a, b, S, B)
#PHI <- sampler(y, mu.0, tau.20, sigma.20, v.0, a, b, S, B)
#plot(mcmc(PHI[, 1:5]))

# drop burnins
#PHI <- PHI[1001:nrow(PHI), ]

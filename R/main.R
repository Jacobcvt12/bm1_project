# Bayesian Methods I Project

# Mixture model of two component univariate mixture of normals

<<<<<<< HEAD
# libraries
library(Rcpp)

# simulated data
y <- 0.3*rnorm(100, 0, 1) + 0.7*rnorm(100, 3, 2)
=======
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
>>>>>>> multisource

# priors
# thetas ~ N(mu, tau^2)
# precisions ~ gamma(v/2, v*sigma^2/2)
# p ~ beta(a, b)

# theta priors
mu.0 <- 0
<<<<<<< HEAD
tau.20 <- 10
=======
tau.20 <- 5
>>>>>>> multisource

# s priors
sigma.20 <- 1
v.0 <- 5

# p priors
<<<<<<< HEAD
a <- 1
b <- 3
=======
a <- 2
b <- 1
>>>>>>> multisource

# load C++ MCMC sampler
#sourceCpp("src/sampler.cpp")

# number of iterations
S <- 1e5
<<<<<<< HEAD

# explore joint posterior
PHI <- sampler(y, mu.0, tau.20, sigma.20, v.0, a, b, S)
=======
B <- 1e4

# explore joint posterior
set.seed(42)
#library(coda)
#PHI <- sampler(y, mu.0, tau.20, sigma.20, v.0, a, b, S, B)[(B+1):(B+S), ]
#plot(mcmc(PHI[, 1:5]))
>>>>>>> multisource

# drop burnins
PHI <- PHI[1001:nrow(PHI), ]

library(MASS)

create.data <- function(params,
                        dir,
                        distribution=rnorm, 
                        output.dims=1,
                        mix.prob=c(0.5, 0.5), 
                        obs=100)
{
    # get number of components
    k <- length(mix.prob)

    # allocate data
    data <- matrix(0, nrow=obs, ncol=output.dims)

    for (i in 1:k) {
        data <- data + mix.prob[k]*distribution(obs, params[k, ])
    }

    write.matrix(data, dir)
}

params <- matrix(c(0, 1,
                   10, 20),
                 byrow=TRUE,
                 nrow=2)

dir <- "data/simulated.dat"

create.data(params, dir)

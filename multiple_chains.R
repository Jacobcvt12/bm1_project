library(doParallel)
#need to detect number of cores to see how many chains you can run simultaneously
cores=detectCores()
#As a rule of thumb, probably want to use cores-2 so you can still do other things
#while the simulation is running
cl=makeCluster(cores-2,type="SOCK")
registerDoParallel(cl)
#number of chains to run in parallel
num_chains=4

chains = foreach(k=1:4,.combine=cbind) %dopar%
         #the argument .combine=cbind will combine the results from each, 
         #the first set of columns will be from chain 1 and so on 
         #going to assume that you will run 4 chains
         {set.seed(1)
          chain1=rnorm(100)
          set.seed(2)
          chain2=rnorm(100)
          set.seed(3)
          chain3=rnorm(100)
          set.seed(4)
          chain4=rnorm(100)
         }  
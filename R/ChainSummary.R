### chain <- 1 chain of MCMC output with columns Mean, Var, K

library(coda)

MEAN_auto_corr<-acf(chain$Mean, plot=FALSE)
MEAN_EffectiveSS<-effectiveSize(chain$Mean)
MEAN_HPD<-HPDinterval(as.mcmc(chain$Mean), prob=0.95)

VAR_auto_corr<-acf(chain$Var, plot=FALSE)
VAR_EffectiveSS<-effectiveSize(chain$Var)
VAR_HPD<-HPDinterval(as.mcmc(chain$Var), prob=0.95)

K_Table<-table(chain$K)

MEAN_auto_corr
MEAN_EffectiveSS
MEAN_HPD

VAR_auto_corr
VAR_EffectiveSS
VAR_HPD

K_table

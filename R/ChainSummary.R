### chain <- 1 chain of MCMC output with columns $P, $Theta1, $Theta2, $s1, $s2, $x1, $x2

library(coda)

P_auto_corr<-acf(chain$P, plot=FALSE)
P_EffectiveSS<-effectiveSize(chain$P)
P_HPD<-HPDinterval(as.mcmc(chain$P), prob=0.95)

THETA1_auto_corr<-acf(chain$Theta1, plot=FALSE)
THETA1_EffectiveSS<-effectiveSize(chain$Theta1)
THETA1_HPD<-HPDinterval(as.mcmc(chain$Theta1), prob=0.95)

S1_auto_corr<-acf(chain$S1, plot=FALSE)
S1_EffectiveSS<-effectiveSize(chain$S1)
S1_HPD<-HPDinterval(as.mcmc(chain$S1), prob=0.95)

THETA2_auto_corr<-acf(chain$Theta2, plot=FALSE)
THETA2_EffectiveSS<-effectiveSize(chain$Theta2)
THETA2_HPD<-HPDinterval(as.mcmc(chain$Theta2), prob=0.95)

S2_auto_corr<-acf(chain$S2, plot=FALSE)
S2_EffectiveSS<-effectiveSize(chain$S2)
S2_HPD<-HPDinterval(as.mcmc(chain$S2), prob=0.95)


pX1<-sum(chain[chain$x1==1]/length(chain$x1)
pX2<-sum(chain[chain$x2==1]/length(chain$x2)


P_auto_corr
P_EffectiveSS
P_HPD

THETA1_auto_corr
THETA1_EffectiveSS
THETA1_HPD

S1_auto_corr
S1_EffectiveSS
S1_HPD<

THETA2_auto_corr
THETA2_EffectiveSS
THETA2_HPD

S2_auto_corr
S2_EffectiveSS
S2_HPD


pX1
pX2

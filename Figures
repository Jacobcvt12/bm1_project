############# Figures

phi1<-PHI[PHI[,106]==1,]
phi2<-PHI[PHI[,106]==2,]

## histogram and ac of theta1 k=1
hist(phi1[,1], xlab=expression(theta[1]), main=expression(paste("Histogram of ", theta[1], " for ", k==1, " model. ", n==34652)), breaks=1000)

theta1k1acf<-acf(phi1[,1], plot=FALSE)
plot(theta1k1acf, main=expression(paste("Autocorrelation plot for ", theta[1],". ", k==1, " model")))

##histogram  and ac of precision1 k=1
hist(phi1[,3], xlab=expression(1/sigma[1]), main=expression(paste("Histogram of ", 1/sigma[1], " for ", k==1, " model. ", n==34652)), breaks=1000)

prec1k1acf<-acf(phi1[,3], plot=FALSE)
plot(prec1k1acf,main=expression(paste("Autocorrelation plot for ", 1/sigma[1],". ", k==1, " model")))

## histogram and ac of theta1 k=2
hist(phi2[,1], xlab=expression(theta[1]), main=expression(paste("Histogram of ", theta[1], " for ", k==2, " model. ", n==65348)), breaks=1000)

theta1k2acf<-acf(phi2[,1], plot=FALSE)
plot(theta1k2acf,main=expression(paste("Autocorrelation plot for ", theta[1],". ", k==2, " model")))

## histogram and ac of precision1 k=2

hist(phi2[,3], xlab=expression(1/sigma[1]), main=expression(paste("Histogram of ", 1/sigma[1], " for ", k==2, " model. ", n==65348)), breaks=1000)

prec1k2acf<-acf(phi2[,3], plot=FALSE)
plot(prec1k2acf,main=expression(paste("Autocorrelation plot for ", 1/sigma[1],". ", k==2, " model")))

## histogram and ac of theta2 k=2
hist(phi2[,2], xlab=expression(theta[2]), main=expression(paste("Histogram of ", theta[2], " for ", k==2, " model. ", n==65348)), breaks=1000)

theta2k2acf<-acf(phi2[,2], plot=FALSE)
plot(theta2k2acf,main=expression(paste("Autocorrelation plot for ", theta[2],". ", k==2, " model")))

## histogram and ac of precision2 k=2

hist(phi2[,4], xlab=expression(1/sigma[2]), main=expression(paste("Histogram of ", 1/sigma[2], " for ", k==2, " model. ", n==65348)), breaks=1000)

prec2k2acf<-acf(phi2[,4], plot=FALSE)
plot(prec2k2acf,main=expression(paste("Autocorrelation plot for ", 1/sigma[2],". ", k==2, " model")))

## histogram and ac of p k=2

hist(phi2[,5], xlab=expression(p), main=expression(paste("Histogram of binomial ", p, " for ", k==2, " model. ", n==65348)), breaks=1000)

pk2acf<-acf(phi2[,5], plot=FALSE)
plot(pk2acf,main=expression(paste("Autocorrelation plot for binomial ", p,". ", k==2, " model")))

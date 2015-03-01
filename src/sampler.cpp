#include <Rcpp.h>
#include <cmath>

// [[Rcpp::export]]
Rcpp::NumericMatrix sampler(Rcpp::NumericVector y, 
                            Rcpp::NumericVector mu0,
                            Rcpp::NumericVector tau20,
                            Rcpp::NumericVector sigma20,
                            Rcpp::NumericVector v0,
                            double a, double b, int S=1000) {
    // allocate PHI matrix for posterior distribution
    Rcpp::NumericMatrix PHI(S, 5 + y.length());

    // set random number generator
    Rcpp::RNGScope scope;

    // simulate first value from prior distribution
    double p = Rcpp::as<double>(Rcpp::rbeta(1, a, b));
    PHI(0, 0) = p;

    double theta1= Rcpp::as<double>(Rcpp::rnorm(1, mu0(0), sqrt(tau20(0))));
    PHI(0, 1) = theta1;

    double s1 = Rcpp::as<double>(Rcpp::rgamma(1, v0(0)/2, 
                                              2/(v0(0)*sigma20(0))));
    PHI(0, 2) = s1;

    double theta2 = Rcpp::as<double>(Rcpp::rnorm(1, mu0(1), sqrt(tau20(1))));
    PHI(0, 3) = theta2;

    double s2 = Rcpp::as<double>(Rcpp::rgamma(1, v0(1)/2, 
                                              2/(v0(1)*sigma20(1))));
    PHI(0, 4) = s1;

    // still need to initialize xi's

    // initialize value for delta
    double delta = 1;
    
    // begin metropolis routine
    for (int s = 1; s < S; ++s) {
        double p = PHI(s-1, 0);
        double p_star = Rcpp::as<double>(Rcpp::rnorm(1, p, delta^2));
    }

    return PHI;
}

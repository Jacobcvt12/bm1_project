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
    double p = R::rbeta(a, b);
    PHI(0, 0) = p;

    double theta1= R::rnorm(mu0(0), sqrt(tau20(0)));
    PHI(0, 1) = theta1;

    double s1 = R::rgamma(v0(0)/2, 2/(v0(0)*sigma20(0)));
    PHI(0, 2) = s1;

    double theta2 = R::rnorm(mu0(1), sqrt(tau20(1)));
    PHI(0, 3) = theta2;

    double s2 = R::rgamma(v0(1)/2, 2/(v0(1)*sigma20(1)));
    PHI(0, 4) = s1;

    // still need to initialize xi's

    // initialize value for delta
    double delta2 = 1;

    // begin metropolis routine
    for (int s = 1; s < S; ++s) {
        // Get last iteration estimation
        p = PHI(s-1, 0);
        theta1 =  PHI(s-1, 1);
        s1 = PHI(s-1, 2);
        theta2 = PHI(s-1, 3);
        s2 = PHI(s-1, 4);

        // compute candidate values
        double p_star = R::rnorm(p, delta2);
        double theta1_star = R::rnorm(theta1, delta2);
        double s1_star = R::rnorm(s1, delta2);
        double theta2_star = R::rnorm(theta2, delta2);
        double s2_star = R::rnorm(s2, delta2);

        // accept or reject candidate
        // likeiihood
        
        // r
    }

    return PHI;
}

#include <Rcpp.h>
#include <cmath>

// [[Rcpp::export]]
Rcpp::NumericMatrix sampler(Rcpp::NumericVector y, 
                            Rcpp::NumericVector mu0,
                            Rcpp::NumericVector tau20,
                            Rcpp::NumericVector sigma20,
                            Rcpp::NumericVector v0,
                            double a, double b, int S=1000) {
    // Initialize PHI matrix for posterior distribution
    //Rcpp::NumericMatrix PHI(1, 5 + y.length());
    Rcpp::NumericMatrix PHI(1, 5);

    // Set random number generator
    Rcpp::RNGScope scope;
    //Rf_PrintValue(Rcpp::rbeta(1, a, b));

    // Simulate first value from prior distribution
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


    return PHI;
}

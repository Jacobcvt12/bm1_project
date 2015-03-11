#include <Rcpp.h>
#include <cmath>
#include <string>
#include <map>
#include <numeric>

#include "prev_vals.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix sampler(Rcpp::NumericVector y, 
                            double mu0, double tau20, 
                            double sigma20, double v0,
                            double a, double b, 
                            int S=1000, int B=1000) {
    // allocate PHI matrix for posterior distribution
    Rcpp::NumericMatrix PHI(S+B, 5 + y.length() + 1);

    // set random number generator
    Rcpp::RNGScope scope;

    // simulate first value from prior distribution
    double p = R::rbeta(a, b);
    //double p = 0.3;
    PHI(0, 0) = p;

    double theta1= R::rnorm(mu0, sqrt(tau20));
    //double theta1 = 0;
    PHI(0, 1) = theta1;

    //double s1 = R::rgamma(v0/2, 2/(v0*sigma20));
    double s1 = 1;
    PHI(0, 2) = s1;

    double theta2 = R::rnorm(mu0, sqrt(tau20));
    //double theta2 = 2;
    PHI(0, 3) = theta2;

    //double s2 = R::rgamma(v0/2, 2/(v0*sigma20));
    double s2 = 1/9;
    PHI(0, 4) = s2;

    // intial x's~binom(0.5)
    Rcpp::NumericVector x_s = Rcpp::rbinom(y.length(), 1, 0.3);
    for (int i = 5; i < 5 + y.length(); ++i) {
        PHI(0, i) = x_s(i-5);
    }

    params worker(2, y);
    worker.initial(2, Rcpp::NumericVector::create(theta1, theta2),
                   Rcpp::NumericVector::create(s1, s2),
                   x_s, p);
    worker.initial(1, Rcpp::NumericVector::create((theta1+theta2)/2, 0),
                   Rcpp::NumericVector::create((s1, s2)/2, 0),
                   Rcpp::rbinom(y.length(), 1, 0), p);

    // begin metropolis routine
    for (int s = 0; s < S + B; ++s) {
        if ((s % 100) == 0) {
            Rcpp::Rcout << s << std::endl;
        }

        // update PHI
        PHI(s, Rcpp::_) = worker.iteration(((s % 100) == 0 & S < B));
    }

    //delta.erase_buffer();

    return PHI;
}


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

    // intial x's~binom(0.5)
    Rcpp::NumericVector x_s = Rcpp::rbinom(y.length(), 1, 0.5);
    for (int i = 5; i < 5 + y.length(); ++i) {
        PHI(0, i) = x_s(i-5);
    }

    // initialize value for delta
    double delta2 = 1;

    // begin metropolis routine
    for (int s = 1; s < S; ++s) {
        // get last iteration estimation
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

        // this generation of x* may be incorrect
        // should it depend on delta?
        Rcpp::NumericVector x_s_star = Rcpp::rbinom(y.length(), 1, 0.5);

        // accept or reject candidate
        double log_r;

        // theta1 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y, theta1_star, 1/sqrt(s1), true)) + 
                 R::dnorm(theta1_star, mu0(0), sqrt(tau20(0)), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y, theta1, 1/sqrt(s1), true)) +
                 R::dnorm(theta1, mu0(0), sqrt(tau20(0)), true));

        if (log(R::runif(0, 1)) < log_r) {
            theta1 = theta1_star;
        }
        
        // s1 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y, theta1, 1/sqrt(s1_star), true)) + 
                 R::dgamma(s1_star, v0(0)/2, 2/(v0(0)*sigma20(0)), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y, theta1, 1/sqrt(s1), true)) +
                 R::dgamma(s1, v0(0)/2, 2/(v0(0)*sigma20(0)), true));

        if (log(R::runif(0, 1)) < log_r) {
            s1 = s1_star;
        }
        
        // theta2 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y, theta2_star, 1/sqrt(s2), true)) + 
                 R::dnorm(theta2_star, mu0(1), sqrt(tau20(1)), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y, theta2, 1/sqrt(s2), true)) +
                 R::dnorm(theta2, mu0(1), sqrt(tau20(1)), true));

        if (log(R::runif(0, 1)) < log_r) {
            theta2 = theta2_star;
        }

        // s2 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y, theta2, 1/sqrt(s2_star), true)) + 
                 R::dgamma(s2_star, v0(0)/2, 2/(v0(0)*sigma20(0)), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y, theta2, 1/sqrt(s2), true)) +
                 R::dgamma(s2, v0(0)/2, 2/(v0(0)*sigma20(0)), true))

        if (log(R::runif(0, 1)) < log_r) {
            s2 = s2_star;
        }
    }

    return PHI;
}

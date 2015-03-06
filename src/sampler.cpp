#include <Rcpp.h>
#include <cmath>
#include <string>
#include <map>
#include <numeric>

#include "track_proposals.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix sampler(Rcpp::NumericVector y, 
                            double mu0, double tau20, 
                            double sigma20, double v0,
                            double a, double b, 
                            int S=1000, int B=1000) {
    // allocate PHI matrix for posterior distribution
    Rcpp::NumericMatrix PHI(S+B, 5 + y.length());

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

    // NB: theta's and s's should be changed to vectors
    // instead of individual variables

    // intial x's~binom(0.5)
    Rcpp::NumericVector x_s = Rcpp::rbinom(y.length(), 1, 0.3);
    for (int i = 5; i < 5 + y.length(); ++i) {
        PHI(0, i) = x_s(i-5);
    }

    // initialize value for delta
    track_acceptance delta;
    delta.track_parameter("p", 0.1);
    delta.track_parameter("theta1", 1);
    delta.track_parameter("s1", 0.2);
    delta.track_parameter("theta2", 1);
    delta.track_parameter("s2", 0.2);
    for (int i = 0; i < x_s.length(); i++) {
        delta.track_parameter("x" + std::to_string (i), 0.3);
    }

    // begin metropolis routine
    for (int s = 1; s < S + B; ++s) {
        if (s == B) {
            // burn in period done
            // following line causes segfault
            // work this out!
            // delta.erase_buffer();
        } else if (s < B) {
            // still in burnin period
            if (s % 100 == 0) {
                delta.modify_deltas();
            }
        }

        // include only y's where corresponding x
        // is appropriate class
        Rcpp::NumericVector y_1 = y[x_s == 1];
        Rcpp::NumericVector y_2 = y[x_s == 0];;

        // propose values
        double p_star = R::rnorm(p, delta["p"]);
        double theta1_star = R::rnorm(theta1, delta["theta1"]);
        double s1_star = R::rnorm(s1, delta["s1"]);
        double theta2_star = R::rnorm(theta2, delta["theta2"]);
        double s2_star = R::rnorm(s2, delta["s2"]);
        Rcpp::NumericVector x_s_star(x_s.length());

        for (int i = 0; i < x_s.length(); i++) {
            if (R::runif(0, 1) < delta["x" + std::to_string (i)]) {
                x_s_star(i) = int (x_s(i) + 1) % 2;
            } else {
                x_s_star(i) = x_s(i);
                delta.accept_reject("x" + std::to_string(i), 0);
            }
        }

        // accept or reject candidate
        double log_r;

        // theta1 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y_1, theta1_star, 1/sqrt(s1), true)) + 
                 R::dnorm(theta1_star, mu0, sqrt(tau20), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_1, theta1, 1/sqrt(s1), true)) +
                 R::dnorm(theta1, mu0, sqrt(tau20), true));

        if (log(R::runif(0, 1)) < log_r) {
            theta1 = theta1_star;
            delta.accept_reject("theta1", 1);
        } else {
            delta.accept_reject("theta1", 0);
        }
        
        // s1 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y_1, theta1, 1/sqrt(s1_star), true)) + 
                 R::dgamma(s1_star, v0/2, 2/(v0*sigma20), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_1, theta1, 1/sqrt(s1), true)) +
                 R::dgamma(s1, v0/2, 2/(v0*sigma20), true));

        if (log(R::runif(0, 1)) < log_r) {
            s1 = s1_star;
            delta.accept_reject("s1", 1);
        } else {
            delta.accept_reject("s1", 0);
        }
        
        // theta2 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y_2, theta2_star, 1/sqrt(s2), true)) + 
                 R::dnorm(theta2_star, mu0, sqrt(tau20), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_2, theta2, 1/sqrt(s2), true)) +
                 R::dnorm(theta2, mu0, sqrt(tau20), true));

        if (log(R::runif(0, 1)) < log_r) {
            theta2 = theta2_star;
            delta.accept_reject("theta2", 1);
        } else {
            delta.accept_reject("theta2", 0);
        }

        // s2 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y_2, theta2, 1/sqrt(s2_star), true)) + 
                 R::dgamma(s2_star, v0/2, 2/(v0*sigma20), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_2, theta2, 1/sqrt(s2), true)) +
                 R::dgamma(s2, v0/2, 2/(v0*sigma20), true));

        if (log(R::runif(0, 1)) < log_r) {
            s2 = s2_star;
            delta.accept_reject("s2", 1);
        } else {
            delta.accept_reject("s2", 0);
        }

        // p log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dbinom(x_s, 1, p_star, true)) +
                 R::dbeta(p_star, a, b, true))
                -
                (Rcpp::sum(Rcpp::dbinom(x_s, 1, p, true)) +
                 R::dbeta(p, a, b, true));

        if (log(R::runif(0, 1)) < log_r) {
            p = p_star;
            delta.accept_reject("p", 1);
        } else {
            delta.accept_reject("p", 0);
        }

        // x log acceptance ratio
        for (int i = 0; i < x_s.length(); i++) {
            if (x_s(i) != x_s_star(i)) {
                if (int (x_s(i)) == 0) {
                    log_r = R::dnorm(y(i), theta1, 1/sqrt(s1), true) +
                            R::dbinom(1, 1, 0.3, true) -
                            R::dnorm(y(i), theta2, 1/sqrt(s2), true) -
                            R::dbinom(0, 1, 0.3, true);
                } else {
                    log_r = R::dnorm(y(i), theta2, 1/sqrt(s2), true) +
                            R::dbinom(0, 1, 0.3, true) -
                            R::dnorm(y(i), theta1, 1/sqrt(s1), true) -
                            R::dbinom(1, 1, 0.3, true);
                } 

                if (log(R::runif(0, 1)) < log_r) {
                    x_s(i) = x_s_star(i);
                    delta.accept_reject("x" + std::to_string(i), 1);
                } else {
                    delta.accept_reject("x" + std::to_string(i), 0);
                }
            } 
        }

        // update PHI
        PHI(s, 0) = p;
        PHI(s, 1) = theta1;
        PHI(s, 2) = s1;
        PHI(s, 3) = theta2;
        PHI(s, 4) = s2;
        
        for (int i = 5; i < 5 + y.length(); i++) {
            PHI(s, i) = x_s(i-5);
        }
    }

    return PHI;
}


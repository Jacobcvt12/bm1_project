#include <Rcpp.h>
#include <cmath>
#include <string>
#include <map>
#include <numeric>
#include <boost/circular_buffer.hpp>

#define ACCEPTTRACK 100

Rcpp::NumericVector y_of_class_n(Rcpp::NumericVector y, 
                                 Rcpp::NumericVector x, 
                                 int k=1)
{
    int n = y.length();
    std::vector<int> which_class;
    which_class.reserve(n);

    for (int i = 0; i < n; i++) {
        if (x(i) == k) {
            which_class.push_back(x(i));
        }
    }
            
    Rcpp::NumericVector y_k(which_class.size());

    for (int i = 0; i < which_class.size(); i++) {
        y_k(i) = y(which_class[i]);
    }

    return y_k;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sampler(Rcpp::NumericVector y, 
                            Rcpp::NumericVector mu0,
                            Rcpp::NumericVector tau20,
                            Rcpp::NumericVector sigma20,
                            Rcpp::NumericVector v0,
                            double a, 
                            double b, 
                            int S=1000,
                            int B=1000) {
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
    std::map<std::string, float> delta;
    delta.insert(std::make_pair("p", 0.1));
    delta.insert(std::make_pair("theta1", 1));
    delta.insert(std::make_pair("s1", 0.1));
    delta.insert(std::make_pair("theta2", 1));
    delta.insert(std::make_pair("s2", 0.1));
                

    // keep running list of last 100 accepted or rejected proposals for each parameter
    std::map<std::string, boost::circular_buffer<int>*> accept_reject;    
    accept_reject.insert(std::make_pair("p", new boost::circular_buffer<int>(ACCEPTTRACK)));
    accept_reject.insert(std::make_pair("theta1", new boost::circular_buffer<int>(ACCEPTTRACK)));
    accept_reject.insert(std::make_pair("s1", new boost::circular_buffer<int>(ACCEPTTRACK)));
    accept_reject.insert(std::make_pair("theta2", new boost::circular_buffer<int>(ACCEPTTRACK)));
    accept_reject.insert(std::make_pair("s2", new boost::circular_buffer<int>(ACCEPTTRACK)));
    
    // This should be updated only through burning stage

    // do i need to manually release the pointers to circular_buffers or will map do that for me?

    // begin metropolis routine
    for (int s = 1; s < S + B; ++s) {
        if (s == B) {
            // burn in done
            // release accept_reject map
            for (auto itr=accept_reject.begin(); itr!=accept_reject.end(); ++itr) {
                delete itr->second;
                accept_reject.erase(itr);
            }
        } else if (s < B) {
            // burn in iterations
            if (s % 100 == 0) {
                for (auto itr=accept_reject.begin(); itr!=accept_reject.end(); ++itr) {
                    // get mean acceptantances
                    double accept_rate = std::accumulate(itr->second.begin(),
                                                         itr->second.end(), 0) / ACCEPTTRACK;

                    // need to do something with this!
                }
            }
        }

        // propose values
        double p_star = R::rnorm(p, delta["p"]);
        double theta1_star = R::rnorm(theta1, delta["theta1"]);
        double s1_star = R::rnorm(s1, delta["s1"]);
        double theta2_star = R::rnorm(theta2, delta["theta2"]);
        double s2_star = R::rnorm(s2, delta["s2"]);

        // this generation of x* may be incorrect
        // how should it depend on delta?
        Rcpp::NumericVector x_s_star = Rcpp::rbinom(y.length(), 1, 0.5);

        // accept or reject candidate
        double log_r;

        // include only y's where corresponding x
        // is appropriate class

        Rcpp::NumericVector y_1 = y_of_class_n(y, x_s, 1);
        Rcpp::NumericVector y_2 = y_of_class_n(y, x_s, 0);

        // theta1 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y_1, theta1_star, 1/sqrt(s1), true)) + 
                 R::dnorm(theta1_star, mu0(0), sqrt(tau20(0)), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_1, theta1, 1/sqrt(s1), true)) +
                 R::dnorm(theta1, mu0(0), sqrt(tau20(0)), true));

        if (log(R::runif(0, 1)) < log_r) {
            theta1 = theta1_star;
            accept_reject["theta1"]->push_back(1);
        } else {
            accept_reject["theta1"]->push_back(0);
        }
        
        // s1 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y_1, theta1, 1/sqrt(s1_star), true)) + 
                 R::dgamma(s1_star, v0(0)/2, 2/(v0(0)*sigma20(0)), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_1, theta1, 1/sqrt(s1), true)) +
                 R::dgamma(s1, v0(0)/2, 2/(v0(0)*sigma20(0)), true));

        if (log(R::runif(0, 1)) < log_r) {
            s1 = s1_star;
        }
        
        // theta2 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y_2, theta2_star, 1/sqrt(s2), true)) + 
                 R::dnorm(theta2_star, mu0(1), sqrt(tau20(1)), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_2, theta2, 1/sqrt(s2), true)) +
                 R::dnorm(theta2, mu0(1), sqrt(tau20(1)), true));

        if (log(R::runif(0, 1)) < log_r) {
            theta2 = theta2_star;
        }

        // s2 log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dnorm(y_2, theta2, 1/sqrt(s2_star), true)) + 
                 R::dgamma(s2_star, v0(0)/2, 2/(v0(0)*sigma20(0)), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_2, theta2, 1/sqrt(s2), true)) +
                 R::dgamma(s2, v0(0)/2, 2/(v0(0)*sigma20(0)), true));

        if (log(R::runif(0, 1)) < log_r) {
            s2 = s2_star;
        }

        // p log acceptance ratio
        log_r = (Rcpp::sum(Rcpp::dbinom(x_s, 1, p_star, true)) +
                 R::dbeta(p_star, a, b, true))
                -
                (Rcpp::sum(Rcpp::dbinom(x_s, 1, p, true)) +
                 R::dbeta(p, a, b, true));

        if (log(R::runif(0, 1)) < log_r) {
            p = p_star;
        }

        // x log acceptance ratio
        for (int i = 0; i < x_s.length(); i++) {
            log_r = (R::dbinom(x_s_star(i), 1, p, true) +
                     R::dbinom(x_s_star(i), 1, 0.5, true))
                    -
                    (R::dbinom(x_s(i), 1, p, true) +
                     R::dbinom(x_s(i), 1, 0.5, true));
            if (log(R::runif(0, 1)) < log_r) {
                x_s(i) = x_s_star(i);
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


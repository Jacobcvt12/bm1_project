#include "prev_vals.h"

params::params(int dimension, Rcpp::NumericVector data) {
    max_dim = dimension;
    curr_dim = dimension;
    y = data;

    delta.track_parameter("p", 0.1);
    delta.track_parameter("theta20", 1);
    delta.track_parameter("s20", 0.2);
    delta.track_parameter("theta21", 1);
    delta.track_parameter("s21", 0.2);
    delta.track_parameter("theta10", 1);
    delta.track_parameter("s10", 0.2);
    for (int i = 0; i < y.length(); i++) {
        delta.track_parameter("x" + std::to_string (i), 0.3);
    }
}

void params::initial(int dimension, 
                     Rcpp::NumericVector theta,
                     Rcpp::NumericVector s,
                     Rcpp::NumericVector x_s,
                     double p)
{
    theta_map.insert(std::make_pair(dimension, theta));
    s_map.insert(std::make_pair(dimension, s));
    x_s_map.insert(std::make_pair(dimension, x_s));
    p_map.insert(std::make_pair(dimension, p));
}

void params::priors(double prior1, 
                    double prior2, 
                    double prior3,
                    double prior4, 
                    double prior5, 
                    double prior6)
{
    mu0=prior1;
    tau20=prior2;
    sigma20=prior3;
    v0=prior4;
    a=prior5;
    b=prior6;
}

void params::update_theta(int update_dim) {
    Rcpp::NumericVector theta_curr = theta_map[update_dim];
    Rcpp::NumericVector s_curr = s_map[update_dim];
    Rcpp::NumericVector x_s_curr = x_s_map[update_dim];
    double p_curr = p_map[update_dim];

    double log_r;

    // loop over dimensions to update each one
    for (int i = 0; i < update_dim; i++) {
        double theta = theta_curr(i);
        double s = s_curr(i);
        Rcpp::NumericVector x_s = x_s_map[update_dim];
        key = "theta" + std::to_string(update_dim) + std::to_string(i);

        // propose theta*
        double theta_star = R::rnorm(theta, delta[key]);

        // subset y by component
        Rcpp::NumericVector y_sub = y[x_s == i];

        //Rcpp::Rcout << "likelihood of theta* " <<
            //Rcpp::sum(Rcpp::dnorm(y_sub, theta_star, 1/sqrt(s), true)) <<
            //std::endl;
        //Rcpp::Rcout << "likelihood of theta " <<
            //Rcpp::sum(Rcpp::dnorm(y_sub, theta, 1/sqrt(s), true)) <<
            //std::endl;
        //Rcpp::Rcout << "likelihood of theta* given priors " <<
            //R::dnorm(theta_star, mu0, sqrt(tau20), true) <<
            //std::endl;
        //Rcpp::Rcout << "likelihood of theta given priors " <<
            //R::dnorm(theta, mu0, sqrt(tau20), true) <<
            //std::endl;

        // check log r
        log_r = (Rcpp::sum(Rcpp::dnorm(y_sub, theta_star, 1/sqrt(s), true)) + 
                 R::dnorm(theta_star, mu0, sqrt(tau20), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_sub, theta, 1/sqrt(s), true)) +
                 R::dnorm(theta, mu0, sqrt(tau20), true));

        if (log(R::runif(0, 1)) < log_r) {
            theta_curr(i) = theta_star;
            delta.accept_reject(key, 1);
        } else {
            delta.accept_reject(key, 0);
        }
    }
    
    theta_map[update_dim] = theta_curr;
}

void params::update_s(int update_dim) {
    Rcpp::NumericVector theta_curr = theta_map[update_dim];
    Rcpp::NumericVector s_curr = s_map[update_dim];
    Rcpp::NumericVector x_s_curr = x_s_map[update_dim];
    double p_curr = p_map[update_dim];

    double log_r;
    double theta, s;

    // loop over dimensions to update each one
    for (int i = 0; i < update_dim; i++) {
        double theta = theta_curr(i);
        double s = s_curr(i);
        Rcpp::NumericVector x_s = x_s_map[update_dim];
        key = "s" + std::to_string(update_dim) + std::to_string(i);

        // propose theta*
        // need to include delta somehow
        double s_star = R::rnorm(s, delta[key]);

        // subset y by component
        Rcpp::NumericVector y_sub = y[x_s == i];

        // check log r
        log_r = (Rcpp::sum(Rcpp::dnorm(y_sub, theta, 1/sqrt(s_star), true)) + 
                 R::dgamma(s_star, v0/2, 2/(v0*sigma20), true))
                -
                (Rcpp::sum(Rcpp::dnorm(y_sub, theta, 1/sqrt(s), true)) +
                 R::dgamma(s, v0/2, 2/(v0*sigma20), true));

        if (log(R::runif(0, 1)) < log_r) {
            s_curr(i) = s_star;
            delta.accept_reject(key, 1);
        } else {
            delta.accept_reject(key, 0);
        }
    }
    
    s_map[update_dim] = s_curr;
}

void params::update_x_s(int update_dim) {
    Rcpp::NumericVector theta_curr = theta_map[update_dim];
    Rcpp::NumericVector s_curr = s_map[update_dim];
    Rcpp::NumericVector x_s_curr = x_s_map[update_dim];
    double p_curr = p_map[update_dim];

    double log_r;
    double theta, s;

    // this is a manual not-so-nice update
    if (update_dim ==2) {
        // propose new x_s
        Rcpp::NumericVector x_s_star(x_s_curr.length());

        for (int i = 0; i < x_s_curr.length(); i++) {
            if (R::runif(0, 1) < delta["x" + std::to_string (i)]) {
                x_s_star(i) = int (x_s_curr(i) + 1) % 2;
            } else {
                x_s_star(i) = x_s_curr(i);
                delta.accept_reject("x" + std::to_string(i), 0);
            }
        }

        // check log r
        for (int i = 0; i < x_s_curr.length(); i++) {
            if (x_s_curr(i) != x_s_star(i)) {
                if (int (x_s_curr(i)) == 0) {
                    log_r = R::dnorm(y(i), theta_curr(1), 
                                     1/sqrt(s_curr(1)), true) +
                            R::dbinom(1, 1, 0.3, true) -
                            R::dnorm(y(i), theta_curr(0), 
                                     1/sqrt(s_curr(0)), true) -
                            R::dbinom(0, 1, 0.3, true);
                } else {
                    log_r = R::dnorm(y(i), theta_curr(1), 
                                     1/sqrt(s_curr(1)), true) +
                            R::dbinom(0, 1, 0.3, true) -
                            R::dnorm(y(i), theta_curr(0), 
                                     1/sqrt(s_curr(0)), true) -
                            R::dbinom(1, 1, 0.3, true);
                } 

                if (log(R::runif(0, 1)) < log_r) {
                    x_s_curr(i) = x_s_star(i);
                    delta.accept_reject("x" + std::to_string(i), 1);
                } else {
                    delta.accept_reject("x" + std::to_string(i), 0);
                }
            } 
        }

        x_s_map[update_dim] = x_s_curr;
    }
}

void params::update_p(int update_dim) {
    Rcpp::NumericVector x_s_curr = x_s_map[update_dim];
    double p_curr = p_map[update_dim];

    double log_r;

    // less sophisticated if then
    if (update_dim == 2) {
        // propose p*
        double p_star = R::rnorm(p_curr, delta["p"]);
        
        log_r = (Rcpp::sum(Rcpp::dbinom(x_s_curr, 1, p_star, true)) +
                 R::dbeta(p_star, a, b, true))
                -
                (Rcpp::sum(Rcpp::dbinom(x_s_curr, 1, p_curr, true)) +
                 R::dbeta(p_curr, a, b, true));

        if (log(R::runif(0, 1)) < log_r) {
            p_curr = p_star;
            delta.accept_reject("p", 1);
        } else {
            delta.accept_reject("p", 0);
        }
    }

    p_map[update_dim] = p_curr;
}

void params::update_k()
{
    // with prob 0.5 propose new k
    if (R::runif(0, 1) < 0.5) {
        int new_dim;
        if (curr_dim == 1) {
            new_dim = 2;
        } else {
            new_dim = 1;
        }

        double likelihood1 = 0;
        double likelihood2 = 0;

        Rcpp::NumericVector theta2 = theta_map[2];
        Rcpp::NumericVector theta1 = theta_map[1];

        Rcpp::NumericVector s2 = s_map[2];
        Rcpp::NumericVector s1 = s_map[1];

        Rcpp::NumericVector x_s2 = x_s_map[2];
        
        double p = p_map[2];

        for (int i = 0; i < y.length(); i++) {
            likelihood1 += R::dnorm(y(i), theta1(0), 
                                    1/sqrt(s1(0)), true);
            if (x_s2(i) == 0) {
                likelihood2 += R::dnorm(y(i), theta2(1), 
                                        1/sqrt(s2(1)), true);
            } else {
                likelihood2 += R::dnorm(y(i), theta2(0), 
                                        1/sqrt(s2(0)), true);
            }
        }

        // downweight likelihood of 1 dimension
        likelihood1 -= 160;

        //Rcpp::Rcout << "likelihood of 1 component " << likelihood1 <<
            //std::endl;
        //Rcpp::Rcout << "likelihood of 2 component " << likelihood2 <<
            //std::endl;
        
        double log_r;

        if (new_dim == 1) {
            log_r = likelihood1 - likelihood2;
        } else {
            log_r = likelihood2 - likelihood1;
        }

        if (log(R::runif(0, 1)) < log_r) {
            curr_dim = new_dim;
        }
    }
}

Rcpp::NumericVector params::iteration(bool update_delta) {
    // initialize row
    Rcpp::NumericVector s(5+y.length()+1);

    // updates
    for (int i = 1; i <= 2; i++) {
        update_theta(i);
        update_s(i);
        update_x_s(i);
        update_p(i);
        update_k();
    }

    s(0) = theta_map[curr_dim](0);
    s(1) = theta_map[curr_dim](1);
    s(2) = s_map[curr_dim](0);
    s(3) = s_map[curr_dim](1);
    s(4) = p_map[curr_dim];

    for (int i = 5; i < y.length(); i++) {
        s(i) = x_s_map[curr_dim](i-5);
    }

    s(5+y.length()) = curr_dim;

    if (update_delta) {
        delta.modify_deltas();
    }

    return s;
}


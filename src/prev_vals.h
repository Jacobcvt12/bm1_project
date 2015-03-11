#ifndef PREV_VALS_H
#define PREV_VALS_H

#include <Rcpp.h>
#include <map>
#include <cmath>

#include "track_proposals.h"

class params {
    // store ys
    Rcpp::NumericVector y;

    // track deltas
    track_acceptance delta;
    std::string key;

    // maps of dimension and current values for given parameter
    std::map<int, Rcpp::NumericVector> theta_map;
    std::map<int, Rcpp::NumericVector> s_map;
    std::map<int, Rcpp::NumericVector> x_s_map;
    std::map<int, double> p_map;

    // information about dimensions
    int max_dim;
    int curr_dim;

    // priors
    double mu0, tau20, sigma20, v0, a, b;

    public:
        params(int dimension, Rcpp::NumericVector data);
        void initial(int dimension, 
                     Rcpp::NumericVector theta,
                     Rcpp::NumericVector s,
                     Rcpp::NumericVector x_s,
                     double p);

        void priors(double prior1, 
                    double prior2, 
                    double prior3,
                    double prior4, 
                    double prior5, 
                    double prior6);

        // change these to not return anything
        void update_theta();
        void update_s();
        void update_x_s();
        void update_p();
        void update_k();
        
        Rcpp::NumericVector iteration(bool update_delta=false);
};

#endif

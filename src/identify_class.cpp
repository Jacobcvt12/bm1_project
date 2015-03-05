#include "identify_class.h"

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

#ifndef IDENTIFY_CLASS_H
#define IDENTIFY_CLASS_H

#include <Rcpp.h>

Rcpp::NumericVector y_of_class_n(Rcpp::NumericVector y, 
                                 Rcpp::NumericVector x, 
                                 int k=1);

#endif

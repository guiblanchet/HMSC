#ifndef sampleNum_h
#define sampleNum_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::vec sampleNum(Rcpp::NumericVector x,
                        int size,
                        bool replace,
                        Rcpp::NumericVector prob);
#endif
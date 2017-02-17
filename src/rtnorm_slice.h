#ifndef rtnorm_slice_h
#define rtnorm_slice_h

#include <RcppArmadillo.h>
#include <Rcpp.h>

typedef Rcpp::NumericMatrix::iterator mat_iterator;
typedef Rcpp::NumericVector::iterator vec_iterator;

double rtnorm_slice(int iter, double mean, double a, double b);

#endif

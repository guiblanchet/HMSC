#ifndef rwishart_h
#define rwishart_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat rwishart(unsigned int df, const arma::mat& S);

#endif

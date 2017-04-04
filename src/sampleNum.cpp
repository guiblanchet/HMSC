#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "sampleNum.h"

using namespace Rcpp ;

arma::vec sampleNum(Rcpp::NumericVector x,
                        int size,
                        bool replace,
                        Rcpp::NumericVector prob) {
  arma::vec ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

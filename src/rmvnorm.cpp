#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rmvnorm.h"

using namespace arma; 
using namespace Rcpp;

//' @title Sample from a multivariate Normal distribution.
//' 
//' @description Produces one or more samples from the specified multivariate normal distribution. Unlike the function \code{\link{mvrnorm}[MASS]}, this function is very simple and does not perform the many checks carried out in the \code{\link{mvrnorm}[MASS]}. For this reason, it is faster but it must be used carefully as it may lead to some unexpected error. This function is only meant to be used internally.
//'
//' @param n Numeric. Number of samples.
//' @param Mean A vector giving the means of the variables sampled. 
//' @param Var A positive-definite symmetric matrix specifying the covariance matrix of the variables.
//'
//' @export
//[[Rcpp::export]]
arma::mat rmvnorm(int n, arma::vec Mean, arma::mat Var) {
   int ncols = Var.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(Mean, 1, n).t() + Y * arma::chol(Var);
}

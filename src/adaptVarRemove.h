#ifndef adaptVarRemove_h
#define adaptVarRemove_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::field<arma::mat> adaptVarRemove(arma::mat& param,
						   arma::mat& varia,
						   arma::mat& shrinkLocal,
						   arma::mat& paramShrinkGlobal,
						   arma::vec& shrinkGlobal,
						   arma::mat& shrink,
						   arma::vec& redund,
						   arma::uvec& notToShrink,
						   double toShrinkSums,
						   double nsp,
						   double nparam,
						   int nVariaLev);

#endif

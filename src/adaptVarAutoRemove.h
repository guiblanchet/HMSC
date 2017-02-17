#ifndef adaptVarAutoRemove_h
#define adaptVarAutoRemove_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::field<arma::mat> adaptVarAutoRemove(arma::mat& param,
										  arma::mat& varia,
										  arma::vec& paramAuto,
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
#ifndef adaptVarAutoAdd_h
#define adaptVarAutoAdd_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::field<arma::mat> adaptVarAutoAdd(arma::mat& param,
									   arma::mat& varia,
									   arma::vec& paramAuto,
									   arma::mat& shrinkLocal,
									   arma::mat& paramShrinkGlobal,
									   arma::vec& shrinkGlobal,
									   arma::mat& shrink,
									   arma::vec& redund,
									   double priorShrinkLocal,
									   double priorShrinkSpeedShape,
									   double priorShrinkSpeedScale,
									   arma::mat& priorParamAutoDist,
									   double nsp,
									   double nparam,
									   int nVariaLev);

#endif
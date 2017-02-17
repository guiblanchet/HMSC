#ifndef updateParamShrinkGlobal_h
#define updateParamShrinkGlobal_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::vec updateParamShrinkGlobal(arma::mat& shrinkLocal,
								  arma::mat& param2,
								  arma::vec& paramShrinkGlobal,
								  arma::vec& shrinkGlobal,
								  double priorShrinkOverallShape,
								  double priorShrinkOverallScale,
								  double priorShrinkSpeedShape,
								  double priorShrinkSpeedScale,
								  double nsp,
								  int nVaria);

#endif

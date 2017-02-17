#ifndef updateMeansParamXBase_h
#define updateMeansParamXBase_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat updateMeansParamXBase(arma::mat& priorMeansParamX, arma::mat& paramX,
									 double nsp, int nparamX);

#endif

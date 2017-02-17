#ifndef updatePrecX_h
#define updatePrecX_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rwishart.h"

arma::mat updatePrecX(arma::mat& meansParamX,
					  arma::mat& priorVarXScaleMat,
					  double priorVarXDf, 
					  arma::mat& paramX, 
					  arma::mat& precX,
					  double nsp);

#endif

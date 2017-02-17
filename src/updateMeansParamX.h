#ifndef updateMeanParamX_h
#define updateMeanParamX_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateMeansParamXBase.h"
#include "updatePrecX.h"

arma::mat updateMeansParamX(arma::mat& priorMeansParamX, 
						   arma::mat& priorVarXScaleMat,
						   double priorVarXDf,
						   arma::mat& paramX, 
						   arma::mat& meansparamX,
						   arma::mat& precX, 
						   double nsp, 
						   int nparamX);

#endif

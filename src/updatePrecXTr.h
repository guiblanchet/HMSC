#ifndef updatePrecXTr_h
#define updatePrecXTr_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rwishart.h"

arma::mat updatePrecXTr(arma::mat& Tr,
						arma::mat& priorVarXScaleMat,
						double priorVarXDf, 
						arma::mat& paramTr, 
						arma::mat& paramX, 
						arma::mat& precX,
						double nsp, 
						int nparamX);

#endif

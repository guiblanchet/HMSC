#ifndef updatePrecXPhylo_h
#define updatePrecXPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rwishart.h"

arma::mat updatePrecXPhylo(arma::mat& meansParamX,
						   arma::mat& paramX,
						   arma::mat& precX,
						   arma::mat& wPhyloInvMat,
						   arma::mat& priorVarXScaleMat,
						   double priorVarXDf,
						   double nsp);

#endif

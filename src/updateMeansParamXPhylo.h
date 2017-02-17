#ifndef updateMeansParamXPhylo_h
#define updateMeansParamXPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rmvnorm.h"

arma::mat updateMeansParamXPhylo(arma::mat& paramX,
								 arma::mat& meansParamX,
								 arma::mat& precX,
								 arma::mat& wPhyloInvMat,
								 arma::mat& priorMeansParamX,
								 arma::mat& priorVarMeansParamX,
								 double nsp);

#endif

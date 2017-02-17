#ifndef updateParamX_h
#define updateParamX_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rmvnorm.h"

arma::mat updateParamX(arma::mat& Ylatent,
					   arma::mat& X, 
					   arma::mat& meansparamX,
					   arma::mat& precX,
					   arma::mat& paramX,
					   double nsp,
					   int nsite,
					   int nparamX);

#endif

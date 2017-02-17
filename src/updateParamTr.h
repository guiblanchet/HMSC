#ifndef updateParamTr_h
#define updateParamTr_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rmvnorm.h"

arma::mat updateParamTr(arma::mat& Tr,
						arma::mat& paramX,
						arma::mat& paramTr,
						arma::mat& precX,
						arma::mat& priorVarTr,
						arma::mat& priorParamTr,
						int nparamX,
						int nTr);

#endif

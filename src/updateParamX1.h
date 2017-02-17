#ifndef updateParamX1_h
#define updateParamX1_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rmvnorm.h"

arma::mat updateParamX1(arma::mat& Ylatent,
						arma::mat& X, 
						arma::mat& meansparamX,
						arma::mat& precX,
						arma::mat& paramX,
						arma::vec& residVar,
						double nsp, 
						int nsite, 
						int nparamX);

#endif

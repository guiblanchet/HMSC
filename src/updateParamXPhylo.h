#ifndef updateParamXPhylo_h
#define updateParamXPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rmvnorm.h"

arma::mat updateParamXPhylo(arma::mat& Ylatent,
							arma::mat& X,
							arma::mat& paramX,
							arma::mat& meansParamX,
							arma::mat& precX,
							arma::vec& residVar,
							arma::mat& wPhyloInvMat,
							double nsp,
							int nsite,
							int nparamX);

#endif

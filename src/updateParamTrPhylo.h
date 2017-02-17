#ifndef updateParamTrPhylo_h
#define updateParamTrPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

arma::mat updateParamTrPhylo(arma::mat& Tr,
							 arma::mat& paramX,
							 arma::mat& paramTr,
							 arma::mat& precX,
							 arma::mat& wPhyloInvMat,
							 arma::mat& priorParamTr,
							 arma::mat& priorVarTr,
							 int nTr,
							 int nparamX);

#endif
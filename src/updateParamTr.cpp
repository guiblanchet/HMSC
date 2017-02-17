#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamTr.h"

using namespace arma;
using namespace Rcpp;

// This function updates paramTr
arma::mat updateParamTr(arma::mat& Tr,
						arma::mat& paramX,
						arma::mat& paramTr,
						arma::mat& precX,
						arma::mat& priorVarTr,
						arma::mat& priorParamTr,
						int nparamX,
						int nTr){

	// Sample from a standard normal distribution to fill paramTrsmpl
//	mat paramTrsmpl=randn(nTr,nparamX);

	//Update paramTr
//	paramTr=paramTr+trans(trans(chol(varTr))*paramTrsmpl*chol(varX));
	
	// Define the variance of the multivariate normal distribution from which paramTr will be sampled
	mat varParamTr = inv(inv(priorVarTr)+kron(precX,Tr*trans(Tr)));
	
	// Define the mean of the multivariate normal distribution from which paramTr will be sampled
	mat meanParamTr = varParamTr*(inv(priorVarTr)*trans(priorParamTr)+vectorise(Tr*paramX*precX));
	
	// Sample from multivariate normal distribution
	paramTr = trans(reshape(rmvnorm(1,meanParamTr,varParamTr),nTr,nparamX));
	
	// Return result
	return paramTr;
}
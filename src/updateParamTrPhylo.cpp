#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateMeansParamXPhylo.h"

using namespace arma;
using namespace Rcpp;

// This function updates paramTr when information on phylogeny is given
arma::mat updateParamTrPhylo(arma::mat& Tr,
							 arma::mat& paramX,
							 arma::mat& paramTr,
							 arma::mat& precX,
							 arma::mat& wPhyloInvMat,
							 arma::mat& priorParamTr,
							 arma::mat& priorVarTr,
							 int nTr,
							 int nparamX){
	
	// Calculate trait weighted wPhyloInv
	mat wPhyloInvTr = Tr*wPhyloInvMat;

	// Calculate variance to sample paramTr
	mat varParamTr = inv(inv(priorVarTr)+kron(precX,wPhyloInvTr*trans(Tr)));

	// Recalculate mean to sample paramTr
	mat meanParamTr = varParamTr*(inv(priorVarTr)*trans(priorParamTr)+vectorise(wPhyloInvTr*paramX*precX));
	
	// Update paramTr
	paramTr = trans(reshape(rmvnorm(1,meanParamTr,varParamTr),nTr,nparamX));
	
	// Return a list of results
	return paramTr;
}
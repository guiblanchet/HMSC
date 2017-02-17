#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateMeansParamXPhylo.h"

using namespace arma;
using namespace Rcpp;

// This function updates the average of the parameters of X (meansParamX)
arma::mat updateMeansParamXPhylo(arma::mat& paramX,
								 arma::mat& meansParamX,
								 arma::mat& precX,
								 arma::mat& wPhyloInvMat,
								 arma::mat& priorMeansParamX,
								 arma::mat& priorVarMeansParamX,
								 double nsp){
	
	// Sum of all vallues in wPhyloInvMat (this is equivalent to 1*wPhyloInvMat*t^t, where "1" = vector of 1s nsp times and t is the transpose of a matrix)
	double wPhyloInvAccu = accu(wPhyloInvMat);
	
	// Calculate variance of meansParamX
	mat varMeansParamX = inv(inv(priorVarMeansParamX)+precX*wPhyloInvAccu);
	
	// Column sums of wPhyloInvMat (this is equivalent to 1*wPhyloInvMat, where "1" = vector of 1s nsp times)
	mat wPhyloInvSum = sum(wPhyloInvMat);
	
	// Recalculate mean to sample meansParamX
	mat meanMeansParamX = varMeansParamX*(inv(priorVarMeansParamX)*priorMeansParamX+trans(wPhyloInvSum*paramX*precX));
	
	// Update meansParamX
	meansParamX = trans(rmvnorm(1,meanMeansParamX,varMeansParamX));
	
	// Return a list of results
	return meansParamX;
}
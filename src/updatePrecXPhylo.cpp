#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updatePrecXPhylo.h"

using namespace arma;
using namespace Rcpp;

// This function updates the precision matrix precX when traits are accounted for
arma::mat updatePrecXPhylo(arma::mat& meansParamX,
						   arma::mat& paramX,
						   arma::mat& precX,
						   arma::mat& wPhyloInvMat,
						   arma::mat& priorVarXScaleMat,
						   double priorVarXDf,
						   double nsp){
		
	// Calculate the amount of information in paramX not explained by traits
	mat residParamX = paramX-meansParamX;
	
	// Calculate the scale matrix to be used for sampling varX from an inverse Wishart distribution
	mat varXScaleMat = inv(priorVarXScaleMat+trans(residParamX)*wPhyloInvMat*residParamX);
	
	// Update precX and varX
	precX = rwishart(priorVarXDf+nsp,varXScaleMat);
	
	// Return results
	return precX;
}

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updatePrecXTr.h"

using namespace arma;
using namespace Rcpp;

// This function updates the precision matrix precX when traits are accounted for
arma::mat updatePrecXTr(arma::mat& Tr,
						arma::mat& priorVarXScaleMat,
						double priorVarXDf, 
						arma::mat& paramTr, 
						arma::mat& paramX, 
						arma::mat& precX,
						double nsp, 
						int nparamX){
		
	// Calculate the amount of information in paramX not explained by traits
	mat residParamX = paramX-trans(paramTr*Tr);
	
	// Calculate the scale matrix to be used for sampling varX from an inverse Wishart distribution
	mat varXScaleMat = inv(priorVarXScaleMat+trans(residParamX)*residParamX);
	
	// Update precX and varX
	precX = rwishart(priorVarXDf+nsp,varXScaleMat);
	
	// Return results
	return precX;
}

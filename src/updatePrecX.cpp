#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updatePrecX.h"

using namespace arma;
using namespace Rcpp;

// This function updates the precision matrix precX

arma::mat updatePrecX(arma::mat& meansParamX,
					  arma::mat& priorVarXScaleMat,
					  double priorVarXDf, 
					  arma::mat& paramX, 
					  arma::mat& precX,
					  double nsp){
		
	// Defining objects
	mat residParamX = paramX-trans(meansParamX*ones(1,nsp));
	
	// Calculate the scale matrix to be used for sampling varX from an inverse Wishart distribution
	mat varXScaleMat = inv(priorVarXScaleMat+trans(residParamX)*residParamX);
	
	// Update precX and varX
	precX = rwishart(priorVarXDf+nsp,varXScaleMat);
	
	// Return results
	return precX;
}
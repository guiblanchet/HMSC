#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateMeansParamX.h"

using namespace arma;
using namespace Rcpp;

// This function updates the average of the parameters of X (meansParamX)
arma::mat updateMeansParamX(arma::mat& priorMeansParamX, 
						   arma::mat& priorVarXScaleMat,
						   double priorVarXDf,
						   arma::mat& paramX, 
						   arma::mat& meansparamX,
						   arma::mat& precX, 
						   double nsp, 
						   int nparamX){

	// Defining objects
	mat varXScaleMat(nparamX,nparamX);
	vec sumParamX(nparamX);
	
	for(int i = 0; i < nparamX; i++){
		sumParamX(i)=sum(paramX.col(i));
	}
	
	mat meansparamXBase(nparamX,1);
	mat meansparamXsmpl(nparamX,1);
	mat paramXprior(nparamX,1);
	
	mat matnsp(nparamX,1);
	matnsp.fill(1/(nsp+1));
	
	// Recalculate meansparamX
	paramXprior = sumParamX+priorMeansParamX;
	meansparamXBase = updateMeansParamXBase(priorMeansParamX,paramX,nsp,nparamX);
	
	// Update meansparamX
	meansparamXsmpl = randn(nparamX);
	meansparamX = meansparamXBase + (sqrt(1/(nsp+1))*(chol(precX.i(),"lower")*meansparamXsmpl));
	
	// Return a list of results
	return meansparamX;
}
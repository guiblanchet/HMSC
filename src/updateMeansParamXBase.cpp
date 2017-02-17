#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateMeansParamXBase.h"

using namespace arma;
using namespace Rcpp;

// This function only calculates meanspParamXBase which is use to updates
//   1 - The average of the parameters of X (meanParamX)
//   2 - The precision matrix precX
//   3 - The variance matrix varX

arma::mat updateMeansParamXBase(arma::mat& priorMeansParamX, arma::mat& paramX,
									 double nsp, int nparamX){

	// Defining objects
	vec sumParamX(nparamX);
	
	for(int i = 0; i < nparamX; i++){
		sumParamX(i)=sum(paramX.col(i));
	}
	
	mat meanspParamXBase(nparamX,1);
	mat paramXprior(nparamX,1);
	
	mat matnsp(nparamX,1);
	matnsp.fill(1/(nsp+1));
	
	// Recalculate meanspParamXBase
	paramXprior = sumParamX+priorMeansParamX;
	meanspParamXBase = paramXprior % matnsp;
	
	// Return results
	return meanspParamXBase;
}
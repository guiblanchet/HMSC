#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "adaptVarRemove.h"

using namespace arma;
using namespace Rcpp;

// This function maybe unnecessary

arma::field<arma::mat> adaptVarRemove(arma::mat& param,
						   arma::mat& varia,
						   arma::mat& shrinkLocal,
						   arma::mat& paramShrinkGlobal,
						   arma::vec& shrinkGlobal,
						   arma::mat& shrink,
						   arma::vec& redund,
						   arma::uvec& notToShrink,
						   double toShrinkSums,
						   double nsp,
						   double nparam,
						   int nVariaLev){
		
	// Defining objects
	int notToShrinkSums = accu(notToShrink); // Number of latent variables not to be removed
			
	// Remove variable(s)
	if(toShrinkSums > 0 && nparam > 2){
		if(toShrinkSums > (nparam - 2)){
			vec nonRedundLatent = regspace(0,1);
		}else{
			uvec nonRedundLatent(notToShrinkSums);
			nonRedundLatent = find(notToShrink);
			nparam = notToShrinkSums;
			varia = varia.cols(nonRedundLatent);
			param = param.cols(nonRedundLatent);;
			shrinkLocal = shrinkLocal.cols(nonRedundLatent);;
			paramShrinkGlobal = paramShrinkGlobal(nonRedundLatent);;
			shrinkGlobal = cumprod(paramShrinkGlobal);
			
			mat shrinkLocalT = trans(shrinkLocal);
			shrink = trans(shrinkLocalT.each_col() % shrinkGlobal);
		}
	}

	//////////////////////////
	// Construct result object
	//////////////////////////
	field<mat> result(7,1);
	
	// New variables
	result(0,0) = varia;
	// New number of parameters
	result(1,0) = nparam;
	// New set parameters
	result(2,0) = param;
	// New shrinkLocal
	result(3,0) = shrinkLocal;
	// New paramShrinkGlobal
	result(4,0) = paramShrinkGlobal;
	// New shrinkGlobal
	result(5,0) = shrinkGlobal;
	// New shrink
	result(6,0) = shrink;
	
	// Return results
	return result;
}
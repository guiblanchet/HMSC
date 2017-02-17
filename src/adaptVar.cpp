#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "adaptVar.h"

using namespace arma;
using namespace Rcpp;

arma::field<arma::mat> adaptVar(arma::mat& param,
					 arma::mat& varia,
					 arma::mat& shrinkLocal,
					 arma::mat& paramShrinkGlobal,
					 arma::vec& shrinkGlobal,
					 arma::mat& shrink,
					 arma::vec& redund,
					 double priorShrinkLocal,
					 double priorShrinkSpeedShape,
					 double priorShrinkSpeedScale,
					 double probAdapt,
					 double nsp,
					 double nparam,
					 int nVariaLev,
					 int count){

	// Defining objects
	umat toShrinkCheckBase = sum(abs(param) < redund(1)); // Boolean
	vec toShrinkCheck = conv_to<vec>::from(toShrinkCheckBase); // Convert umat -> mat

	vec toShrinkProp = toShrinkCheck/nsp;
	uvec toShrink = toShrinkProp >= redund(0);
	uvec notToShrink = toShrinkProp < redund(0);
	double toShrinkSums = accu(toShrink);
	
	// Define result object
	field<mat> result(7,1);
	
	// Define a random value sampled from a uniform distribution
	vec randVal = randu(1);
	if(randVal(0) <= probAdapt){
		if(count > 20 && toShrinkSums == 0 && all(toShrinkProp < 0.995)){
			// Add variable
			result = adaptVarAdd(param, varia, shrinkLocal, paramShrinkGlobal, shrinkGlobal, shrink, redund, priorShrinkLocal, priorShrinkSpeedShape, priorShrinkSpeedScale, nsp, nparam, nVariaLev);
		}else{
			// Remove variable(s)
			result = adaptVarRemove(param, varia, shrinkLocal, paramShrinkGlobal, shrinkGlobal, shrink, redund, notToShrink, toShrinkSums, nsp, nparam, nVariaLev);
		}
	}else{
		// Present old results if no adaptation is carried out
		result(0,0) = varia;
		result(1,0) = nparam;
		result(2,0) = param;
		result(3,0) = shrinkLocal;
		result(4,0) = paramShrinkGlobal;
		result(5,0) = shrinkGlobal;
		result(6,0) = shrink;
	}
	
	// Return results
	return result;
}
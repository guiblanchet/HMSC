#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamShrinkGlobal.h"

using namespace arma ;
using namespace Rcpp ;

// This function updates shrinkGlobal, paramShrinkGlobal, and shrink

arma::vec updateParamShrinkGlobal(arma::mat& shrinkLocal,
								  arma::mat& param2,
								  arma::vec& paramShrinkGlobal,
								  arma::vec& shrinkGlobal,
								  double priorShrinkOverallShape,
								  double priorShrinkOverallScale,
								  double priorShrinkSpeedShape,
								  double priorShrinkSpeedScale,
								  double nsp,
								  int nVaria){

	// Define various objects
	mat paramw = shrinkLocal%param2;
	mat colSumParamw = trans(sum(paramw,0)); // The "0" is to say that the cumulative sum is carried out by columns
	double shapeParamShrinkGlobal;
	double scaleParamShrinkGlobal;

	// Update overall shrink level
	shapeParamShrinkGlobal = priorShrinkOverallShape+0.5*nsp*nVaria;
	scaleParamShrinkGlobal = priorShrinkOverallScale+0.5*(1/paramShrinkGlobal(0))*sum(shrinkGlobal%colSumParamw);
	paramShrinkGlobal(0)=randg(1,distr_param(shapeParamShrinkGlobal,1/scaleParamShrinkGlobal))(0);

	// Update shrinkage speed
	for (int i = 1; i < nVaria; i++) {
		shapeParamShrinkGlobal = priorShrinkSpeedShape+0.5*nsp*(nVaria-i+2);
		scaleParamShrinkGlobal = 1/(priorShrinkSpeedScale+0.5*(1/paramShrinkGlobal(i))*sum(shrinkGlobal.tail(nVaria-i)%colSumParamw.tail_rows(nVaria-i)));
		paramShrinkGlobal(i) = randg(1,distr_param(shapeParamShrinkGlobal,scaleParamShrinkGlobal))(0);
		shrinkGlobal = cumprod(paramShrinkGlobal);
	}

	return paramShrinkGlobal;
}

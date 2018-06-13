#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "adaptVarAdd.h"

using namespace arma;
using namespace Rcpp;

arma::field<arma::mat> adaptVarAdd(arma::mat& param,
						arma::mat& varia,
						arma::mat& shrinkLocal,
						arma::mat& paramShrinkGlobal,
						arma::vec& shrinkGlobal,
						arma::mat& shrink,
						arma::vec& redund,
						double priorShrinkLocal,
						double priorShrinkSpeedShape,
						double priorShrinkSpeedScale,
						double nsp,
						double nparam,
						int nVariaLev){

	// Add a variable
	nparam = nparam+1;
	varia = join_rows(varia,randn(nVariaLev));
	param = join_rows(param, mat(nsp,1,fill::zeros));
	shrinkLocal = join_rows(shrinkLocal,randg(nsp,distr_param(priorShrinkLocal/2,2/priorShrinkLocal)));
	paramShrinkGlobal = join_cols(paramShrinkGlobal,randg(1,distr_param(priorShrinkSpeedShape,priorShrinkSpeedScale)));
	shrinkGlobal = cumprod(paramShrinkGlobal);

	mat shrinkLocalT = trans(shrinkLocal);
	shrink = trans(shrinkLocalT.each_col() % shrinkGlobal);

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

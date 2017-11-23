#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateShrinkLocal.h"

using namespace arma ;
using namespace Rcpp ;

arma::mat updateShrinkLocal(arma::mat& shrinkLocal,
							double priorShrinkLocal,
							arma::vec& shrinkGlobal,
							arma::mat& param2,
							double nsp,
							int nparam){

	// Construct transpose of param2
	mat param2T = trans(param2);

	// Update ShrinkLocal
	mat rateShrinkLocal = 1/trans(0.5*(param2T.each_col() % shrinkGlobal)+priorShrinkLocal/2);

	for (int i = 0; i < nsp; i++) {
		for (int j = 0; j < nparam; j++) {
			shrinkLocal(i,j) = randg(1,distr_param(priorShrinkLocal/2+0.5,rateShrinkLocal(i,j)))(0);
		}
	}
	return shrinkLocal;
}

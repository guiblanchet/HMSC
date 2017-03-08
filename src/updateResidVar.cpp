#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateResidVar.h"

using namespace arma ;
using namespace Rcpp ;

// This function updates the inverse residual variance parameter of the model
arma::vec updateResidVar(arma::mat& Yresid,
						arma::vec& residVar,
						double priorResidVarShape,
						double priorResidVarScale,
						double nsp,
						double nsite){

	// Update residVar
	double shapeParam = priorResidVarShape+nsite/2;
	double scaleParam;

	for (int i = 0; i < nsp; i++) {
		scaleParam = (priorResidVarScale+sum(square(Yresid.col(i)))/2);
		residVar(i) = as_scalar(randg(1,distr_param(shapeParam,1/scaleParam)));
	}

	return residVar;
}

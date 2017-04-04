#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateVarDist.h"

using namespace arma ;
using namespace Rcpp ;

arma::vec updateVarDist(arma::mat& Yresid,
						arma::vec& varDist,
						arma::vec& priorVarDistShape,
						arma::vec& priorVarDistScale,
						double nsite,
						double nsp){

	// square all elements of Yresid
	mat Yresid2 = square(Yresid);

	// Update varDist
	for (int i = 0; i < nsp; i++) {
		double shape = priorVarDistShape(i)+nsite/2;
		double scale = 1/(priorVarDistScale(i)+accu(Yresid2.col(i))/2);
		varDist(i) = as_scalar(1/randg(1,distr_param(shape,scale)));
	}

	return varDist;
}

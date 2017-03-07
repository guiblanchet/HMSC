#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamX.h"

using namespace arma ;
using namespace Rcpp ;

arma::mat updateParamX(arma::mat& Ylatent,
					   arma::mat& X,
					   arma::mat& meansparamX,
					   arma::mat& precX,
					   arma::mat& paramX,
					   double nsp,
					   int nsite,
					   int nparamX){

	// Defining objects
	mat varXEst(nparamX,nparamX);
	vec meansparamXEst(nparamX);
	mat Xt = trans(X);

	// Update paramX
	varXEst = inv(precX+Xt*X);
	for (int i = 0; i < nsp; i++) {
		meansparamXEst = varXEst*(precX*trans(meansparamX.row(i))+Xt*Ylatent.col(i));
		paramX.row(i) = rmvnorm(1,meansparamXEst,varXEst);
	}

	return paramX;
}

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamX1.h"

using namespace arma ;
using namespace Rcpp ;

arma::mat updateParamX1(arma::mat& Ylatent,
						arma::mat& X, 
						arma::mat& meansparamX,
						arma::mat& precX,
						arma::mat& paramX,
						arma::vec& residVar,
						double nsp, 
						int nsite, 
						int nparamX){
	
	// Defining objects
	mat varXEst(nparamX,nparamX);
	vec meansparamXEst(nparamX);
	mat Xt = trans(X);
	mat XtX = Xt*X;
	
	// Update paramX
	for (int i = 0; i < nsp; i++) {
		varXEst = inv(precX+XtX*residVar(i));
		meansparamXEst = varXEst*(precX*trans(meansparamX.row(i))+(Xt*Ylatent.col(i))*residVar(i));
		paramX.row(i) = rmvnorm(1,meansparamXEst,varXEst);
	}
	return paramX;
}
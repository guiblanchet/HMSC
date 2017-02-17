#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamXPhylo.h"

using namespace arma ;
using namespace Rcpp ;

arma::mat updateParamXPhylo(arma::mat& Ylatent,
							arma::mat& X,
							arma::mat& paramX,
							arma::mat& meansParamX,
							arma::mat& precX,
							arma::vec& residVar,
							arma::mat& wPhyloInvMat,
							double nsp,
							int nsite,
							int nparamX){
	
	// Defining objects
	mat varXEst(nparamX,nparamX);
	vec meansparamXEst(nparamX);
	mat XtX = trans(X)*X;
	mat residVarDiag = diagmat(residVar);
	
	// Update paramX
	varXEst = inv(kron(precX,wPhyloInvMat)+kron(XtX,residVarDiag));
	meansparamXEst = varXEst*(vectorise(wPhyloInvMat*meansParamX*precX)+vectorise(residVarDiag*trans(Ylatent)*X));
	
	paramX = reshape(rmvnorm(1,meansparamXEst,varXEst),nsp,nparamX);
	
	return paramX;
}
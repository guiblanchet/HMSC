#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamLatentTheoryMultiSp.h"

using namespace arma ;
using namespace Rcpp ;

arma::field<arma::mat> updateParamLatentTheoryMultiSp(arma::mat& Ylatent,
										 arma::umat& Random,
										 arma::vec& residVar,
										 arma::field<arma::mat>& paramLatent,
										 arma::field<arma::mat>& latent,
										 arma::field<arma::mat>& shrink,
										 int nRandom,
										 arma::vec& nLatent,
										 double nsp,
										 int nsite,
									 	 arma::cube& diagMat){


	// Define various objects
	mat Yresid = Ylatent;

	for (int j = 0; j < nRandom; j++) {
		mat shrinkLatent(nLatent(j),nLatent(j));
		mat cholShrinkLatent(nLatent(j),nLatent(j));
		mat meanparamLatentSp(nLatent(j),1);

		// Transpose paramLatent for the analysis (quicker and solve's a dimension problem)
		mat paramLatentT = trans(paramLatent(j,0));

		// Update paramLatent for each species
		for (int k = 0; k < nsp; k++) {
			mat latentCross = trans(diagMat.slice(k)*latent(j,0).rows(Random.col(j)))*diagMat.slice(k)*latent(j,0).rows(Random.col(j));
			shrinkLatent = diagmat(shrink(j,0).row(k))+residVar(k)*latentCross;
			mat latentSp = residVar(k)*trans(diagMat.slice(k)*latent(j,0).rows(Random.col(j)))*Yresid.col(k);
			cholShrinkLatent = chol(shrinkLatent,"lower");
			meanparamLatentSp = solve(trimatu(trans(cholShrinkLatent)),solve(trimatl(cholShrinkLatent),latentSp));
			mat structLatentSp = solve(trimatu(trans(cholShrinkLatent)),randn(nLatent(j)));
			paramLatentT.col(k) = structLatentSp+meanparamLatentSp;
		}

		// Save paramLatentT into paramLatent object
		paramLatent(j,0) = trans(paramLatentT);
	}
	return paramLatent;

}

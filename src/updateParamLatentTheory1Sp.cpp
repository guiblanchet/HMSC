#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamLatentTheory1Sp.h"

using namespace arma ;
using namespace Rcpp ;

arma::field<arma::mat> updateParamLatentTheory1Sp(arma::mat& Ylatent,
										 arma::umat& Random,
										 arma::vec& residVar,
										 arma::field<arma::mat>& paramLatent,
										 arma::field<arma::mat>& latent,
										 arma::field<arma::mat>& shrink,
										 int nRandom,
										 arma::vec& nLatent,
										 double nsp,
										 int nsite,
									 	 arma::mat& diagMat){


	// Define various objects
	mat Yresid(nsite,nsp);

	// Update paramLatent
	for (int j = 0; j < nRandom; j++) {
		Yresid = Ylatent;

		for (int k = 0; k < nRandom; k++) {
			if(!(j == k)){
				Yresid = Yresid - diagMat*latent(k,0).rows(Random.col(k))*trans(paramLatent(k,0));
			}
		}

		mat shrinkLatent(nLatent(j),nLatent(j));
		mat cholShrinkLatent(nLatent(j),nLatent(j));
		mat meanparamLatentSp(nLatent(j),1);

		// Transpose paramLatent for the analysis (quicker and solve's a dimension problem)
		mat paramLatentT = trans(paramLatent(j,0));

		// Update paramLatent for each species
		mat latentCross = trans(diagMat*latent(j,0).rows(Random.col(j)))*diagMat*latent(j,0).rows(Random.col(j));
		for (int k = 0; k < nsp; k++) {
			shrinkLatent = diagmat(shrink(j,0).row(k))+residVar(k)*latentCross;
			mat latentSp = residVar(k)*trans(diagMat*latent(j,0).rows(Random.col(j)))*Yresid.col(k);
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

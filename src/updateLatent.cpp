#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateLatent.h"

using namespace arma ;
using namespace Rcpp ;

arma::field<arma::mat> updateLatent(arma::mat& Ylatent,
									arma::umat& Random,
									arma::vec residVar,
									arma::field<arma::mat>& paramLatent,
									arma::field<arma::mat>& latent,
									int nRandom,
									arma::vec& nRandomLev,
									arma::vec& nLatent,
									double nsp,
									int nsite){

	// Define various objects
	mat Yresid(nsite,nsp);

	// Update latent
	for (int j = 0; j < nRandom; j++) {
		Yresid = Ylatent;

		for (int k = 0; k < nRandom; k++) {
			if(!(j == k)){
				Yresid = Yresid - latent(k,0).rows(Random.col(k))*trans(paramLatent(k,0));
			}
		}

		// Define various objects
		mat precLatent(nLatent(j),nLatent(j));
		mat Q(nLatent(j),nLatent(j));
		mat RprecLatent(nLatent(j),nLatent(j));
		mat invRprecLatent(nLatent(j),nLatent(j));
		mat varLatent(nLatent(j),nLatent(j));
		mat meansLatent(nsite,nLatent(j));
		double wRandomLev;

		mat diagResidVar(nsp,nsp);
		diagResidVar = diagmat(residVar);

		// If the random effect is at the sampling unit level (faster!)
		if(nRandomLev(j)==nsite){
			// Residual weighted paramLatent
			mat wparamLatent = diagResidVar*paramLatent(j,0);

			// Precision matrix for the latent variables
			precLatent = eye(nLatent(j),nLatent(j))+trans(paramLatent(j,0))*wparamLatent;

			// Extract the R matrix from a QR decomposition of the Cholesky decomposition of the precision matrix
			qr(Q,RprecLatent,chol(precLatent));

			// Inverse the R matrix
			invRprecLatent = inv(RprecLatent);

			// Calculate the variance-covariance matrix to sample the of the inverted R matrix
			varLatent = invRprecLatent*trans(invRprecLatent);

			// Calculate the means for the new latent variables
			meansLatent = Yresid*wparamLatent*varLatent;

			// Calculate the new latent variables
			latent(j,0) = meansLatent+randn(nsite,nLatent(j))*trans(invRprecLatent);

		// For any other random effect
		}else{

			uvec RandomLev = unique(Random.col(j));

			// Residual weighted paramLatent
			mat wparamLatent = diagResidVar*paramLatent(j,0);

			// paramLatent squared weighted by inverse residual variance (residVar is always inverted)
			mat paramLatent2 = trans(paramLatent(j,0))*wparamLatent;

			for(int k = 0; k < nRandomLev(j); k++){
				uvec lev = Random.col(j)==RandomLev(k);

				// Weight to add to consider the importance of each levels
				wRandomLev = sum(lev);

				// Precision matrix for the latent variables
				precLatent = eye(nLatent(j),nLatent(j))+paramLatent2*wRandomLev;

				// Extract the R matrix from a QR decomposition of the Cholesky decomposition
				qr(Q,RprecLatent,chol(precLatent));

				// Inverse the R matrix
				invRprecLatent = inv(RprecLatent);

				// Calculate the variance-covariance matrix to sample the of the inverted R matrix
				varLatent = invRprecLatent*trans(invRprecLatent);

				// Calculate the means for the new latent variables
				meansLatent = sum(Yresid.rows(find(lev)))*wparamLatent*varLatent;

				// Calculate the new latent variables
				latent(j,0).row(k) = meansLatent+(randn(1,nLatent(j))*trans(invRprecLatent));
			}
		}
	}

	return latent;
}

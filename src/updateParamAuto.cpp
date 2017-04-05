#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamAuto.h"

using namespace arma;
using namespace Rcpp;

arma::field<arma::vec> updateParamAuto(arma::field<arma::mat>& latentAuto,
									   arma::field<arma::vec>& paramAuto,
									   int npriorParamAuto,
									   arma::mat wAutoDet,
									   arma::field< arma::cube >& wAutoInv,
									   arma::mat& priorParamAutoWeight,
									   Rcpp::NumericMatrix& priorParamAutoDist,
									   arma::vec& nLatentAuto,
									   int nAuto){

	// Define basic objects
	double quadLikeVal;

	for(int h = 0; h < nAuto; h++){
		// Define some other basic objects
		mat latentAutoMat = latentAuto(h,0);
		cube wAutoInvCube = wAutoInv(h,0);
		mat quadLike(npriorParamAuto,nLatentAuto(h));

		for(int i = 0; i < npriorParamAuto; i++){
			// Calculate the quadratic part of a multivariate normal log-likelihood that will be used to sample paramAuto
			quadLike.row(i) = sum(square(wAutoInvCube.slice(i)*latentAutoMat));
		}

		vec paramAutoTmp(nLatentAuto(h));
		// Construct an object to save likelihood results
		Rcpp::NumericVector like(npriorParamAuto); // Rcpp code style

		for(int i = 0; i < nLatentAuto(h); i++){
			for(int j = 0; j < npriorParamAuto; j++){
				// Extract the right value
				quadLikeVal = as_scalar(quadLike(j,i));

				// Calculate the first part of the multivariate normal log-likelihood
				like(j) = log(priorParamAutoWeight(j))-0.5*wAutoDet(j,h)-0.5*quadLikeVal;
			}

			// Recentre log-likelihood
			like = like-max(like); // Rcpp code style

			// Convert log-likelihood to likelihood
			like = exp(like);

			// Normalize the likelihood
			like = like/sum(like);

			// Sample paramAuto
			paramAutoTmp(i) = sampleNum(priorParamAutoDist,1,true,like)(0);
		}

		// Save object
		paramAuto(h,0) = paramAutoTmp;

	}
	// Return result
	return paramAuto;
}

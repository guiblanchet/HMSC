#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamPhylo.h"

using namespace arma;
using namespace Rcpp;

// This function updates paramPhylo when there are traits
double updateParamPhylo(arma::mat& paramX,
						arma::mat& meansParamX,
						double paramPhylo,
						arma::mat& precX,
						arma::vec& wPhyloDet,
						arma::cube& wPhyloInv,
						int nparamX,
						int nparamPhylo,
						int nsp,
						arma::mat& priorParamPhyloWeight,
						Rcpp::NumericVector& priorParamPhyloGrid){

	// Defining objects
	double likeDetPrecX;
	double quadLike;
	Rcpp::NumericVector like(nparamPhylo); // Rcpp code style
	
	// Remove the influence of traits on paramX
	vec paramXresid = vectorise(paramX - meansParamX);

	// Calculate the probability for a predefined set of paramPhylo
	for (int i = 0; i < nparamPhylo; i++){
		// Calculate the quadratic part of a multivariate normal log-likelihood that will be used to sample paramPhylo
		quadLike = as_scalar(trans(paramXresid)*kron(precX,wPhyloInv.slice(i))*paramXresid);
		
		// Calculate the first part of the multivariate normal log-likelihood
		likeDetPrecX = -nsp*log(det(precX))+nparamX*wPhyloDet(i);
		
		// Calculate multivariate normal log-likelihood
		like(i) = log(priorParamPhyloWeight(i,0))-0.5*likeDetPrecX-0.5*quadLike; 
	}
	
	// Recentre log-likelihood
	like = like-max(like); // Rcpp code style
	
	// Convert log-likelihood to likelihood
	like = exp(like);
	
	// Normalize the likelihood
	like = like/sum(like);
	
	// Sample paramPhylo
	paramPhylo = sampleNum(priorParamPhyloGrid,1,true,like)(0);
	
	// Return result
	return paramPhylo;
}
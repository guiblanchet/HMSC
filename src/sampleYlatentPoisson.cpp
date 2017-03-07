#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentPoisson.h"

using namespace arma ;
using namespace Rcpp ;
//' @title Sample the model response matrix
//'
//' @description Sample the model response matrix after the Poisson link function was applied. This function is meant to be used internally.
//'
//' @param Y The species community matrix to be modelled.
//' @param Ylatent Model site by species community matrix after the link function is applied.
//' @param EstModel Estimated model for the site by species community matrix.
//' @param residVar Vector of parameters with as many entries as there are species. Each values in this vector describes how much variation is not explained by the model. This vector should only contains 1s, but it was included in this function to deal with potential situations that may arise (in other words, in case I forgot a special case).
//' @param nsp Numeric. Number of species in \code{Y}.
//' @param nsite Numeric. Number of site sampled.
//'
//' @export
// [[Rcpp::export]]
arma::mat sampleYlatentPoisson(arma::mat& Y,
							  arma::mat& Ylatent,
							  arma::mat& EstModel,
							  arma::vec residVar,
							  double nsp,
							  int nsite){

	// Basic objects
	mat Yresid(nsite,nsp);
	vec residSd = sqrt(residVar);

	// Define upper and lower boundary of the truncated normal distribution
	mat low = log(Y) - EstModel;
	mat high = log(Y+1) - EstModel;

	// Sample from a truncated normal distribution to calculate the residual of Ylatent
	for (int i = 0; i < nsite; i++) {
		for (int j = 0; j < nsp; j++) {
			Yresid(i,j) = rtnorm(0, residSd(j), low(i,j), high(i,j));
		}
	}

	// Recalculate Ylatent
	Ylatent = Yresid+EstModel;

	// Return Ylatent
	return Ylatent;
}

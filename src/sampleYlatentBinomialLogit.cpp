#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentBinomialLogit.h"

using namespace arma ;
using namespace Rcpp ;
//' @title Sample the model response matrix
//'
//' @description Sample the model response matrix after the Poisson link function was applied. This function is meant to be used internally.
//'
//' @param Y The species community matrix to be modelled.
//' @param Ylatent Model site by species community matrix after the link function is applied.
//' @param EstModel Estimated model for the site by species community matrix.
//' @param ncount Numeric. Maximum number of individual for a species.
//' @param nsp Numeric. Number of species in \code{Y}.
//' @param nsite Numeric. Number of site sampled.
//'
//' @export
// [[Rcpp::export]]
arma::mat sampleYlatentBinomialLogit(arma::mat& Y,
							  arma::mat& Ylatent,
							  arma::mat& EstModel,
							  double ncount,
							  double nsp,
							  int nsite){

	// Basic objects
	mat Yresid(nsite,nsp);

	// Define upper and lower boundary of the truncated normal distribution for all but the largest value
	mat pLow = Y/ncount;
	mat pHigh = (Y+1)/ncount;

	// Define upper and lower boundary of the truncated normal distribution for
	uvec maxY = find(Y==ncount);
	if(maxY.n_elem > 0){
		pLow.elem(maxY) = (ncount-1)/ncount;
		pHigh.elem(maxY) = ncount;
	}

	mat low = log(pLow/(1-pLow)) - EstModel;
	mat high = log(pHigh/(1-pHigh)) - EstModel;

	// Find all the NAs and fill low and high associated to NAs with -Inf and +Inf
	uvec YNALoc = find_nonfinite(Y);
	vec nasPos(ncount(YNALoc));
	nasPos.fill(datum::inf);

	vec nasNeg(ncount(YNALoc));
	nasNeg.fill(-datum::inf);

	low.elem(YNALoc) = nasNeg;
	high.elem(YNALoc) = nasPos;

	// Sample from a truncated normal distribution to calculate the residual of Ylatent
	for (int i = 0; i < nsite; i++) {
		for (int j = 0; j < nsp; j++) {
			Yresid(i,j) = rtnorm(0, residSd(j), low(i,j), high(i,j));
			if(std::abs(Yresid(i,j)) == datum::inf){
				if(std::abs(high(i,j)) == datum::inf){
					Yresid(i,j) = low(i,j);
				}
				if(std::abs(low(i,j)) == datum::inf){
					Yresid(i,j) = high(i,j);
				}else{
					Yresid(i,j) = high(i,j);
				}
			}
		}
	}

	// Recalculate Ylatent
	Ylatent = Yresid+EstModel;

	// Return Ylatent
	return Ylatent;
}

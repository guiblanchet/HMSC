#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentProbit.h"

using namespace arma ;
using namespace Rcpp ;

//' @title Sample the model response matrix
//'
//' @description Sample the model response matrix after the Probit link function was applied. This function is meant to be used internally.
//'
//' @param Y0Loc A vector defining the locations of all 0s in the community matrix (\code{Y}).
//' @param Y1Loc A vector defining the locations of all 1s in the community matrix (\code{Y}).
//' @param Ylatent Model site by species community matrix after the link function is applied.
//' @param EstModel Estimated model for the site by species community matrix.
//' @param residVar Vector of parameters with as many entries as there are species. Each values in this vector describes how much variation is not explained by the model. This vector should only contains 1s, but it was included in this function to deal with potential situations that may arise (in other words, in case I forgot a special case).
//' @param nsp Numeric. Number of species in \code{Y}.
//' @param nsite Numeric. Number of site sampled.
//'
//' @export
// [[Rcpp::export]]
arma::mat sampleYlatentProbit(arma::uvec& Y0Loc,
							  arma::uvec& Y1Loc,
							  arma::mat& Ylatent,
							  arma::mat& EstModel,
							  arma::vec residVar,
							  double nsp,
							  int nsite){

	// Define basic objects
	mat Yresid(nsite,nsp);

	// Calculate residual standard deviation
	vec residSd = sqrt(residVar);
	// Negative Estimation model
	mat EstModelNeg = -EstModel;

	// Define upper and lower boundary of the truncated normal distribution
	mat low(nsite,nsp);
	low.fill(-datum::inf);

	mat high(nsite,nsp);
	high.fill(datum::inf);

	// Fill low and high with the right part of EstModelNeg
	low.elem(Y1Loc) = EstModelNeg.elem(Y1Loc);
	high.elem(Y0Loc) = EstModelNeg.elem(Y0Loc);

	// Sample from a truncated normal distribution to calculate the residual of Ylatent
	for (int i = 0; i < nsite; i++) {
		for (int j = 0; j < nsp; j++) {
			Yresid(i,j) = rtnorm(0, residSd(j), low(i,j), high(i,j));
		}
	}
	// Find the values in Yresid that are equal to infinite
	uvec YresidPosInf = find(Yresid==datum::inf);

	// Replace the values equal to infinite by the value in low+0.5
	Yresid.elem(YresidPosInf) = low.elem(YresidPosInf)+0.5;

	// Find the values in Yresid that are equal to negative infinite
	uvec YresidNegInf = find(Yresid==-datum::inf);

	// Replace the values equal to negative infinite by the value in high-0.5
	Yresid.elem(YresidNegInf) = high.elem(YresidNegInf)-0.5;

	// Recalculate Ylatent
	Ylatent = EstModel+Yresid;

	// Return Ylatent
	return Ylatent;
}

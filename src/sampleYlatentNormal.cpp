#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentNormal.h"

using namespace arma ;
using namespace Rcpp ;

//' @title Sample the model response matrix
//'
//' @description Sample the model response matrix after the iddentity link function was applied. This function is meant to be used internally.
//'
//' @param Y The species community matrix to be modelled.
//' @param EstModel Estimated model for the site by species community matrix.
//' @param residVar Vector of parameters with as many entries as there are species. Each values in this vector describes how much variation is not explained by the model. This vector should only contains 1s, but it was included in this function to deal with potential situations that may arise (in other words, in case I forgot a special case).
//' @param nsp Numeric. Number of species in \code{Y}.
//' @param nsite Numeric. Number of site sampled.
//'
//' @details
//'
//' This function is meant to be used only when there are NAs in the normally distributed response variables.
//'
//' @export
// [[Rcpp::export]]
arma::mat sampleYlatentNormal(arma::mat& Y,
							  arma::mat& EstModel,
							  arma::vec& residVar,
							  double nsp,
							  int nsite){

	// Basic object
	uvec YNALoc = find_nonfinite(Y);
	mat Ylatent = Y;

	// Reconstruct Ylatent
	mat YlatentTmp = EstModel+randn(nsite,nsp)*diagmat(sqrt(residVar)); // Note : sigma*x+mu where sigma = standard deviation, x is sampled from N(0,1), mu = mean

	// Replace NAs in Ylatent by values in YlatentTmp
	Ylatent.elem(YNALoc) = YlatentTmp.elem(YNALoc);
	// Correct for positive extreme values
	uvec YlatentExtPos = find(Ylatent > 20);

	if(YlatentExtPos.n_elem > 0){
		vec fillVecPos(YlatentExtPos.n_elem);
		fillVecPos.fill(20);
		Ylatent.elem(YlatentExtPos) = fillVecPos;
	}

	// Correct for negative extreme values
	uvec YlatentExtNeg = find(Ylatent < -20);

	if(YlatentExtNeg.n_elem > 0){
		vec fillVecNeg(YlatentExtNeg.n_elem);
		fillVecNeg.fill(-20);
		Ylatent.elem(YlatentExtNeg) = fillVecNeg;
	}

	// Return Ylatent
	return Ylatent;
}

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentPoisson.h"

using namespace arma ;
using namespace Rcpp ;
////// Needs to be tested....
////// Needs to be tested....
////// Needs to be tested....
////// Needs to be tested....
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
arma::mat sampleYlatentPoisson(arma::mat& Y,
							  arma::mat& Ylatent,
							  arma::mat& EstModel,
							  arma::vec residVar,
							  double nsp,
							  int nsite){

//=====================
// From the matlab code
//=====================
	// Define basic objects
//	mat Yresid(nsite,nsp);

	// Basis constant
//	double log1000 = log(1000);

	// Calculate Yresid
//	mat Ycenter = Ylatent - log1000;

//	Yresid = abs(Ycenter);
	//=================================================================
	// Variance of the normal distribution from which we sample Ycenter
	//=================================================================
	// Mean of the sampler
//	mat meanSampler = 500 / Yresid % tanh(Yresid/2);

	// Variance of the sampler
//	mat varSampler1 = (sinh(Yresid)-Yresid) / square(cosh(Yresid/2)); // This converge to 2 as values in Yresid get larger

//	varSampler1.replace(datum::nan,2);

//	mat varSampler = 250/pow(Yresid,3) % varSampler1; // 250 = 1000/4

	// sample from a Normal distribution
//	mat sampleNormVar = abs(randn(nsite,nsp) % sqrt(varSampler) + meanSampler);

	// Calcuate YcenterVar
//	mat YcenterVar = 1/(sampleNormVar + 1);

	//=============================================================
	// Mean of the normal distribution from which we sample Ycenter
	//=============================================================
//	mat YcenterMean = YcenterVar % (Y - 500 + (EstModel-log1000));

	//===============
	// Sample Ycenter
	//===============
//	Ycenter = randn(nsite,nsp) % sqrt(YcenterVar) + YcenterMean;

	//====================
	// Recalculate Ylatent
	//====================
//	Ylatent = Ycenter + log1000;
	// Return Ylatent
//	return Ylatent;
//}

//======================
// From the copula paper
//======================
// Define basic objects
	mat Yresid(nsite,nsp);
	vec residSd = sqrt(residVar);

	mat low = log(Y) - EstModel;
	mat high = log(Y+1) - EstModel;

	for (int i = 0; i < nsite; i++) {
		for (int j = 0; j < nsp; j++) {
			Yresid(i,j) = rtnorm(0, residSd(j), low(i,j), high(i,j));
		}
	}

	Ylatent = Yresid+EstModel;

	// Return Ylatent
	return Ylatent;
}

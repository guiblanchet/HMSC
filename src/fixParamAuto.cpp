#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "fixParamAuto.h"

using namespace arma ;
using namespace Rcpp ;

// Calculate autocorrelated parameters that do not change in the Gibbs sampler
arma::field<arma::cube> fixParamAuto(arma::mat& Auto,
									 Rcpp::NumericMatrix& priorParamAutoDist,
									 double nsp,
									 arma::vec& nAutoLev,
									 int npriorParamAuto,
									 int i){

	// Define various objects
	mat eyeSp(nsp,nsp,fill::eye);

	cube wAutoDet(1,1,npriorParamAuto);

	// Autocorrelation weighted by its parameters for each parameter considered
	// Calculate Euclidean distance among pairs of samples
	mat dist = euclDist(Auto);

	// Initiate wAutoInv as a cube
	cube wAutoInv1(nAutoLev(i),nAutoLev(i),npriorParamAuto);

	for (int j = 0; j < npriorParamAuto; j++) {
		if(priorParamAutoDist(j,i)< 0.00001){
			// If the distance considered is small assume that there is no weight
			// The determinant is 1
			wAutoDet.slice(j) = 1;
			// The inverse of the weight matrix is eye()
			wAutoInv1.slice(j) = eye(nAutoLev(i),nAutoLev(i));
		}else{
			// If the distance considered is larger than 0.00001
			mat wAuto = exp(-dist/priorParamAutoDist(j,0));

			// Calculate the determinant
			wAutoDet.slice(j) = 2*sum(log(diagvec(chol(wAuto))));
			// Calculate the inverse of the weight matrix
			wAutoInv1.slice(j) = inv(wAuto);
		}
	}

	//////////////////////////
	// Construct result object
	//////////////////////////
	field<cube> result(2,1);

	// Determinant of wPhylo
	result(0,0) = wAutoDet;

	// Determinant of wPhylo
	result(1,0) = wAutoInv1;

	// Return results
	return result;
}

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "fixParamPhylo.h"

using namespace arma ;
using namespace Rcpp ;

// Calculate Phylo parameters that do not change in the Gibbs sampler
arma::field<arma::cube> fixParamPhylo(arma::mat& Phylo,
									  arma::mat& iPhylo,
									  double paramPhylo,
									  Rcpp::NumericVector& priorParamPhyloGrid,
									  double nsp,
									  int nparamPhylo){

	// Define various objects
	mat eyeSp(nsp,nsp,fill::eye);
	
	mat wPhylo(nsp,nsp);
	mat wPhylo1(nsp,nsp);
	cube wPhyloDet(1,1,nparamPhylo);
	cube wPhyloInv(nsp,nsp,nparamPhylo);
	
	// Phylogeny weighted by its parameters for each grid parameter considered
	for (int i = 0; i < nparamPhylo; i++) {
		// Calculate the first part of the weighted phylogeny
		if(priorParamPhyloGrid(i) >= 0){
			wPhylo1 = priorParamPhyloGrid(i) * Phylo;
		}else{
			wPhylo1 = (-priorParamPhyloGrid(i)) * iPhylo;
		}
		// Weighted phylogeny parameters
		wPhylo = wPhylo1 + (1-std::abs(priorParamPhyloGrid(i)))*eyeSp;
		wPhyloDet.slice(i) = 2*sum(log(diagvec(chol(wPhylo))));
		wPhyloInv.slice(i) = inv(wPhylo);
	}
	
	//////////////////////////
	// Construct result object
	//////////////////////////
	field<cube> result(2,1);
	
	// Determinant of wPhylo
	result(0,0) = wPhyloDet;
	
	// Determinant of wPhylo
	result(1,0) = wPhyloInv;
	
	// Return results
	return result;
}
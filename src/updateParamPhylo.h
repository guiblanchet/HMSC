#ifndef updateParamPhylo_h
#define updateParamPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleNum.h"

// This function updates paramPhylo when traits are considered
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
						Rcpp::NumericVector& priorParamPhyloGrid);

#endif

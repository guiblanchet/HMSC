#ifndef fixParamPhylo_h
#define fixParamPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "fixParamPhylo.h"

arma::field<arma::cube> fixParamPhylo(arma::mat& Phylo,
									  arma::mat& iPhylo,
									  double paramPhylo,
									  Rcpp::NumericVector& priorParamPhyloGrid,
									  double nsp,
									  int nparamPhylo);

#endif
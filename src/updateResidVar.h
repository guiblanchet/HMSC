#ifndef updateResidVar_h
#define updateResidVar_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::vec updateResidVar(arma::mat& Yresid,
						 arma::vec& residVar,
						 double priorResidVarShape,
						 double priorResidVarScale,
						 double nsp,
						 double nsite);

#endif

#ifndef updateVarDist_h
#define updateVarDist_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::vec updateVarDist(arma::mat& Yresid,
						arma::vec& varDist,
						arma::vec& priorVarDistShape,
						arma::vec& priorVarDistScale,
						double nsite,
						double nsp);

#endif
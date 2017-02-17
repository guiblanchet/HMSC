#ifndef adaptVarAdd_h
#define adaptVarAdd_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::field<arma::mat> adaptVarAdd(arma::mat& param,
						arma::mat& varia,
						arma::mat& shrinkLocal,
						arma::mat& paramShrinkGlobal,
						arma::vec& shrinkGlobal,
						arma::mat& shrink,
						arma::vec& redund,
						double priorShrinkLocal,
						double priorShrinkSpeedShape,
						double priorShrinkSpeedScale,
						double nsp,
						double nparam,
						int nVariaLev);

#endif
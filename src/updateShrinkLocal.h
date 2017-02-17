#ifndef updateShrinkLocal_h
#define updateShrinkLocal_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat updateShrinkLocal(arma::mat& shrinkLocal,
							double priorShrinkLocal,
							arma::vec& shrinkGlobal,
							arma::mat& param2,
							double nsp,
							int nparamX);

#endif

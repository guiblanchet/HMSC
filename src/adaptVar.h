#ifndef adaptVar_h
#define adaptVar_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "adaptVarAdd.h"
#include "adaptVarRemove.h"

arma::field<arma::mat> adaptVar(arma::mat& param,
					 arma::mat& varia,
					 arma::mat& shrinkLocal,
					 arma::mat& paramShrinkGlobal,
					 arma::vec& shrinkGlobal,
					 arma::mat& shrink,
					 arma::vec& redund,
					 double priorShrinkLocal,
					 double priorShrinkSpeedShape,
					 double priorShrinkSpeedScale,
					 double probAdapt,
					 double nsp,
					 double nparam,
					 int nVariaLev,
					 int count);

#endif

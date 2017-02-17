#ifndef adaptVarAuto_h
#define adaptVarAuto_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "adaptVarAutoAdd.h"
#include "adaptVarAutoRemove.h"

arma::field<arma::mat> adaptVarAuto(arma::mat& param,
									arma::mat& varia,
									arma::vec& paramAuto,
									arma::mat& shrinkLocal,
									arma::mat& paramShrinkGlobal,
									arma::vec& shrinkGlobal,
									arma::mat& shrink,
									arma::vec& redund,
									double priorShrinkLocal,
									double priorShrinkSpeedShape,
									double priorShrinkSpeedScale,
									arma::mat& priorParamAutoDist,
									double probAdapt,
									double nsp,
									double nparam,
									int nVariaLev,
									int count);

#endif
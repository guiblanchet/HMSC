#ifndef mcmcNormalXTrLatent_h
#define mcmcNormalXTrLatent_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamX1.h"
#include "updatePrecXTr.h"
#include "updateParamTr.h"
#include "updateParamLatent.h"
#include "updateLatent.h"
#include "updateShrinkLocal.h"
#include "updateParamShrinkGlobal.h"
#include "adaptVar.h"
#include "updateResidVar.h"

RcppExport SEXP mcmcNormalXTrLatent(arma::mat& Ylatent,
									arma::mat& X, 
									arma::mat& Tr, 
									arma::umat& Random, 
									arma::mat& paramX,
									arma::mat& paramTr,
									arma::mat& precX,
									arma::vec& residVar,
									arma::mat& priorParamTr,
									arma::mat& priorVarTr,
									arma::mat& priorVarXScaleMat,
									double priorVarXDf,
									double priorResidVarScale,
									double priorResidVarShape,
									double priorShrinkLocal,
									double priorShrinkOverallShape,
									double priorShrinkOverallScale,
									double priorShrinkSpeedShape,
									double priorShrinkSpeedScale,
									arma::vec& adapt,
									arma::vec& redund,
									int nRandom,
									arma::vec& nRandomLev,
									arma::vec& nLatent,
									double nsp,
									int nsite,
									int nparamX,
									int nTr,
									int niter,
									int nburn,
									int thin,
									int verbose);
							
#endif

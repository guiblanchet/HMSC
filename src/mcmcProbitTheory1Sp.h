#ifndef mcmcProbitTheory1Sp_h
#define mcmcProbitTheory1Sp_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentProbit.h"
#include "updateParamLatentTheory1Sp.h"
#include "updateParamAuto.h"
#include "updateLatentAutoTheory1Sp.h"
#include "updateShrinkLocal.h"
#include "updateParamShrinkGlobal.h"
#include "adaptVarAuto.h"
#include "updateParamX.h"
#include "updatePrecX.h"
#include "updateMeansParamX.h"
#include "fixParamAuto.h"

RcppExport SEXP mcmcProbitTheory1Sp(arma::mat& Y,
								arma::mat& Ylatent,
								arma::mat& X,
								arma::field< arma::mat >& Auto,
								arma::umat& RandomAuto,
								arma::mat& paramX,
								arma::mat& meansParamX,
								arma::mat& precX,
								arma::vec& residVar,
								arma::field< arma::vec >& paramAuto,
								arma::field< arma::mat >& latentAuto,
								arma::field< arma::mat >& paramLatentAuto,
								arma::field< arma::mat >& shrinkLocalAuto,
								arma::field< arma::vec >& paramShrinkGlobalAuto,
								arma::mat& priorMeansParamX,
								arma::mat& priorVarXScaleMat,
								double priorVarXDf,
								double priorResidVarScale,
								double priorResidVarShape,
								arma::mat& priorParamAutoWeight,
								Rcpp::NumericMatrix& priorParamAutoDist,
								double priorShrinkLocal,
								double priorShrinkOverallShape,
								double priorShrinkOverallScale,
								double priorShrinkSpeedShape,
								double priorShrinkSpeedScale,
								arma::vec& adapt,
								arma::vec& redund,
								int nAuto,
								arma::vec& nAutoLev,
								arma::vec& nLatentAuto,
								double nsp,
								int nsite,
								int nparamX,
								int npriorParamAuto,
								int niter,
								int nburn,
								int thin,
								int verbose,
								arma::mat& diagMat);

#endif

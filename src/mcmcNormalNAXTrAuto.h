#ifndef mcmcNormalNAXTrPhylo_h
#define mcmcNormalNAXTrPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentNormal.h"
#include "updateParamX.h"
#include "updateParamTr.h"
#include "updateParamLatent.h"
#include "updateParamAuto.h"
#include "updatePrecXTr.h"
#include "updateLatentAuto.h"
#include "updateShrinkLocal.h"
#include "updateParamShrinkGlobal.h"
#include "adaptVarAuto.h"
#include "fixParamAuto.h"
#include "updateResidVar.h"

RcppExport SEXP mcmcNormalNAXTrAuto(arma::mat& Y,
								  arma::mat& Ylatent,
								  arma::mat& X,
								  arma::mat& Tr,
								  arma::field< arma::mat >& Auto,
								  arma::umat& RandomAuto,
								  arma::mat& paramX,
								  arma::mat& paramTr,
								  arma::mat& precX,
								  arma::vec& residVar,
								  arma::field< arma::vec >& paramAuto,
								  arma::field< arma::mat >& latentAuto,
								  arma::field< arma::mat >& paramLatentAuto,
								  arma::field< arma::mat >& shrinkLocalAuto,
								  arma::field< arma::vec >& paramShrinkGlobalAuto,
								  arma::mat& priorParamTr,
								  arma::mat& priorVarTr,
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
								  int nTr,
								  int npriorParamAuto,
								  int niter,
								  int nburn,
								  int thin,
								  int verbose);

#endif

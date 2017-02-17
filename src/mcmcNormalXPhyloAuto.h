#ifndef mcmcNormalXPhyloAuto_h
#define mcmcNormalXPhyloAuto_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamXPhylo.h"
#include "updateMeansParamXPhylo.h"
#include "updatePrecXPhylo.h"
#include "fixParamPhylo.h"
#include "updateParamPhylo.h"
#include "updateParamLatent.h"
#include "updateParamAuto.h"
#include "updateLatentAuto.h"
#include "updateShrinkLocal.h"
#include "updateParamShrinkGlobal.h"
#include "adaptVarAuto.h"
#include "fixParamAuto.h"
#include "updateResidVar.h"

RcppExport SEXP mcmcNormalXPhyloAuto(arma::mat& Ylatent,
									 arma::mat& X, 
									 arma::mat& Phylo,
									 arma::mat& iPhylo,
									 arma::field< arma::mat >& Auto,
									 arma::umat& RandomAuto,
									 arma::mat& paramX,
									 arma::mat& meansParamX,
									 arma::mat& precX,
									 double paramPhylo,
									 arma::vec& residVar,
									 arma::field< arma::vec >& paramAuto,
									 arma::field< arma::mat >& latentAuto,
									 arma::field< arma::mat >& paramLatentAuto,
									 arma::field< arma::mat >& shrinkLocalAuto,
									 arma::field< arma::vec >& paramShrinkGlobalAuto,
									 arma::mat& priorMeansParamX,
									 arma::mat& priorVarMeansParamX,
									 arma::mat& priorVarXScaleMat,
									 double priorVarXDf,
									 double priorResidVarScale,
									 double priorResidVarShape,
									 arma::mat& priorParamPhyloWeight,
									 Rcpp::NumericVector& priorParamPhyloGrid,
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
									 int nparamPhylo,
									 int npriorParamAuto,
									 int niter,
									 int nburn,
									 int thin,
									 int verbose);

#endif

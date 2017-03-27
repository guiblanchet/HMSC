#ifndef mcmcOverPoissonXTrPhylo_h
#define mcmcOverPoissonXTrPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentPoisson.h"
#include "updateParamXPhylo.h"
#include "updatePrecXPhylo.h"
#include "updateParamTrPhylo.h"
#include "updateParamPhylo.h"
#include "fixParamPhylo.h"
#include "updateResidVar.h"

RcppExport SEXP mcmcOverPoissonXTrPhylo(arma::mat& Y,
								   arma::mat& Ylatent,
								   arma::mat& X,
								   arma::mat& Tr,
								   arma::mat& Phylo,
								   arma::mat& iPhylo,
								   arma::mat& paramX,
								   arma::mat& paramTr,
								   double paramPhylo,
								   arma::mat& precX,
								   arma::vec& residVar,
								   arma::mat& priorParamTr,
								   arma::mat& priorVarTr,
								   arma::mat& priorVarXScaleMat,
								   double priorVarXDf,
								   double priorResidVarScale,
								   double priorResidVarShape,
								   arma::mat& priorParamPhyloWeight,
								   Rcpp::NumericVector& priorParamPhyloGrid,
								   double nsp,
								   int nsite,
								   int nparamX,
								   int nTr,
								   int nparamPhylo,
								   int niter,
								   int nburn,
								   int thin,
								   int verbose);

#endif

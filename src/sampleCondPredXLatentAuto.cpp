#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPredXLatentAuto.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//[[Rcpp::export]]
arma::cube sampleCondPredXLatentAuto(arma::mat& Y,
					 arma::mat& X,
					 arma::field< arma::mat >& Auto,
					 arma::umat& Random,
					 arma::umat& RandomAuto,
					 arma::cube& paramX,
					 arma::field< arma::mat >& latent,
					 arma::field< arma::mat >& paramLatent,
					 arma::field< arma::mat >& latentAuto,
					 arma::field< arma::mat >& paramLatentAuto,
					 arma::field<arma::vec>& paramAuto,
					 arma::mat residVar,
					 Rcpp::NumericMatrix& priorParamAutoDist,
					 int nsite,
					 double nsp,
					 int nRandom,
					 arma::vec& nRandomLev,
					 int nAuto,
					 arma::vec& nAutoLev,
					 int npriorParamAuto,
					 int niter,
					 int nsample,
				 	 std::string family) {

	// Define basic objects to store results
	mat YlatentSample(nsite, nsp);
	cube Ylatent(nsite, nsp, nsample);
	mat EstModel(nsite, nsp);
	EstModel.zeros();

	mat Yresid(nsite,nsp);

	mat residVarT = trans(residVar);

	// Define field object to store one iterations of latent variables and their associated parameters
	field<mat> latent1iter(nRandom,1);
	field<mat> paramLatent1iter(nRandom,1);
	vec nLatent(nRandom);

	field<mat> latentAuto1iter(nAuto,1);
	field<mat> paramLatentAuto1iter(nAuto,1);
	vec nLatentAuto(nAuto);

	//Define wAutoInv object
	field<cube> wAutoInv(nAuto,1);
	field<cube> paramAutoNoGibb(2,1);
	for (int i = 0; i < nAuto; i++) {
		paramAutoNoGibb = fixParamAuto(Auto(i,0),priorParamAutoDist,nsp,nAutoLev,npriorParamAuto,i);
		wAutoInv(i,0) = paramAutoNoGibb(1,0);
	}

	// Convert NumericMatrix to arma::mat
	mat priorParamAutoDistArma = as<arma::mat>(priorParamAutoDist);

	// Reorganize latent
	field<mat> latentOrg(niter, nRandom);
	field<mat> paramLatentOrg(niter, nRandom);
	int counter = 0;
	for(int j = 0; j < nRandom ; j++){
		for(int i = 0; i < niter ; i++){
			latentOrg(i,j) = latent(counter,0);
			paramLatentOrg(i,j) = paramLatent(counter,0);
			counter++;
		}
	}

	// Reorganize latentAuto
	field<mat> latentAutoOrg(niter, nAuto);
	field<mat> paramLatentAutoOrg(niter, nAuto);
	counter = 0;
	for(int j = 0; j < nAuto ; j++){
		for(int i = 0; i < niter ; i++){
			latentAutoOrg(i,j) = latentAuto(counter,0);
			paramLatentAutoOrg(i,j) = paramLatentAuto(counter,0);
			counter++;
		}
	}

	////////////////////////////////
	// Sample conditional prediction
	////////////////////////////////
	for(int i = 0; i < niter ; i++){
		// Calculate the model estimation
		EstModel = X * trans(paramX.slice(i));

		for(int j = 0; j < nRandom ; j++){
			mat latentMat = latentOrg(i,j);
			nLatent(j) = latentMat.n_cols;
			mat paramLatentMat = paramLatentOrg(i,j);
			EstModel = EstModel + latentMat.rows(Random.col(j))*trans(paramLatentMat);
		}

		for(int j = 0; j < nAuto ; j++){
			mat latentAutoMat = latentAutoOrg(i,j);
			nLatentAuto(j) = latentAutoMat.n_cols;
			mat paramLatentAutoMat = paramLatentAutoOrg(i,j);
			EstModel = EstModel + latentAutoMat.rows(RandomAuto.col(j))*trans(paramLatentAutoMat);
		}

		/////////////////
		// Sample Ylatent
		/////////////////
		for (int j = 0; j < nsample; j++) {
			if(family == "probit"){
				uvec Y0Loc = find(Y==0);
				uvec Y1Loc = find(Y==1);
				uvec YNALoc = find_nonfinite(Y);

				mat YlatentIni = zeros(nsite,nsp);

				YlatentSample = sampleYlatentProbit(Y0Loc, Y1Loc, YNALoc, YlatentIni, EstModel, residVarT.col(i), nsp, nsite);
			}

			if(family == "gaussian"){
				mat repResidVar(nsite,nsp);
				repResidVar = repmat(residVarT.col(i),nsite,1);
				YlatentSample = randn(nsite,nsp)%repResidVar+EstModel;
			}

			if(family == "poisson" | family == "overPoisson"){
				mat YlatentIni = zeros(nsite,nsp);
				YlatentSample = sampleYlatentPoisson(Y, YlatentIni, EstModel, residVarT.col(i), nsp, nsite);
			}

			////////////////
			// Update latent
			////////////////
			for(int j = 0; j < nRandom ; j++){
				latent1iter(j,0) = latentOrg(i,j);
				paramLatent1iter(j,0) = paramLatentOrg(i,j);
			}

			for(int j = 0; j < nAuto ; j++){
				latentAuto1iter(j,0) = latentAutoOrg(i,j);
				paramLatentAuto1iter(j,0) = paramLatentAutoOrg(i,j);
			}

			// Remove influence of X variables
			Yresid = YlatentSample - X * trans(paramX.slice(i));

			// Remove influence of latentAuto
			for(int j = 0; j < nAuto ; j++){
				Yresid = Yresid - latentAuto1iter(j,0) * trans(paramLatentAuto1iter(j,0));
			}

			latent1iter = updateLatent(Yresid,Random,residVarT.col(i),paramLatent1iter,latent1iter,nRandom,nRandomLev,nLatent,nsp,nsite);

			for(int j = 0; j < nRandom ; j++){
				latentOrg(i,j) = latent1iter(j,0);
			}

			////////////////////
			// Update latentAuto
			////////////////////

			// Remove influence of X variables
			Yresid = YlatentSample - X * trans(paramX.slice(i));

			// Remove influence of latentAuto
			for(int j = 0; j < nRandom ; j++){
				Yresid = Yresid - latent1iter(j,0) * trans(paramLatent1iter(j,0));
			}

			latentAuto1iter = updateLatentAuto(Yresid, RandomAuto, residVarT.col(i), paramAuto, wAutoInv, paramLatentAuto1iter, latentAuto1iter, priorParamAutoDistArma, nAuto, nAutoLev, nLatentAuto, nsp, nsite);

			for(int j = 0; j < nAuto ; j++){
				latentAutoOrg(i,j) = latentAuto1iter(j,0);
			}
		}
		// Save results
		Ylatent.slice(i) = YlatentSample;
	}

	// return result object
	return Ylatent;
}

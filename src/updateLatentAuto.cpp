#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateLatentAuto.h"

using namespace arma ;
using namespace Rcpp ;

arma::field<arma::mat> updateLatentAuto(arma::mat& Yresid,
										arma::umat& RandomAuto,
										arma::vec& residVar,
										arma::field<arma::vec>& paramAuto,
										arma::field<arma::cube>& wAutoInv,
										arma::field<arma::mat>& paramLatentAuto,
										arma::field<arma::mat>& latentAuto,
										arma::mat& priorParamAutoDistArma,
										int nAuto,
										arma::vec& nAutoLev,
										arma::vec& nLatentAuto,
										double nsp,
										int nsite){

	// Initiate a few basic object
	int beginRow;
	int beginCol;
	int endRow;
	int endCol;

	// Update latentAuto
	for (int j = 0; j < nAuto; j++) {

		// Remove the effect of other autocorrelated latent variables, if any
		for (int k = 0; k < nAuto; k++) {
			if(j != k){
				Yresid = Yresid - latentAuto(k,0).rows(RandomAuto.col(k))*trans(paramLatentAuto(k,0));
			}
		}

		// Define various objects
		field<arma::mat> wAutoInvSel(1,nLatentAuto(j));
		vec seqVec = regspace<vec>(1,nAutoLev(j));

		// Extract the cube of wAutoInv associated to the factor of Auto considered
		cube wAutoInvFactor = wAutoInv(j,0);

		// Extract the matrix used to make the estimation in wAutoInv and convert it to a vector
		for(int k = 0; k < nLatentAuto(j); k++){
			// Find which paramAuto in priorParamAutoDistArma is the one used
			uvec paramAutoPointer = find(priorParamAutoDistArma.col(j)==as_scalar(paramAuto(j,0)(k)));
			wAutoInvSel(0,k) = wAutoInvFactor.slice(paramAutoPointer(0));
		}

		// Construct a matrix of 0
		mat wAutoInvDiag = zeros(nLatentAuto(j)*nAutoLev(j),nLatentAuto(j)*nAutoLev(j));

		// Fill wAutoInvDiag on the diagonal with the matrix selected in wAutoInvSel
		for(int k = 1; k <= nLatentAuto(j); k++){
			beginRow = (k-1)*nAutoLev(j);
			beginCol = (k-1)*nAutoLev(j);
			endRow = (k-1)*nAutoLev(j)+nAutoLev(j)-1;
			endCol = (k-1)*nAutoLev(j)+nAutoLev(j)-1;

			wAutoInvDiag.submat(beginRow,beginCol,endRow,endCol) = wAutoInvSel(0,k-1);
		}

		// Construct a diagonal matrix with residVar
		mat diagResidVar(nsp,nsp);
		diagResidVar = diagmat(residVar);
		///////////////////////////////////////////////////////////////
		// If the random effect is at the sampling unit level (faster!)
		///////////////////////////////////////////////////////////////
		if(nAutoLev(j)==nsite){
			// Residual weighted paramLatentAuto;
			mat wparamLatentAuto = diagResidVar*paramLatentAuto(j,0);

			// Precision matrix for the latentAuto variables
			mat precLatent = wAutoInvDiag + kron(trans(paramLatentAuto(j,0))*wparamLatentAuto,eye(nsite,nsite));

			// Calculate the variance-covariance matrix for the latentAuto variables
			mat varLatent = inv_sympd(precLatent);

			// Calculate the means for the new latentAuto variables
			mat meansLatent = solve(precLatent,vectorise(Yresid * wparamLatentAuto));

			// Calculate the new latentAuto variables
			latentAuto(j,0) = reshape(rmvnorm(1, meansLatent,varLatent), nAutoLev(j), nLatentAuto(j));
		//////////////////////////////
		// For any other random effect
		//////////////////////////////
		}else{
			// Residual weighted paramLatentAuto;
			mat wparamLatentAuto = diagResidVar * paramLatentAuto(j,0);

			// Construct a matrix defining how to weight each sample
			mat diagMat = eye(nAutoLev(j),nAutoLev(j));
			uvec lev = vectorise(RandomAuto.col(j));
			mat diagMatLev = diagMat.cols(lev);

			// Precision matrix for the latentAuto variables
			mat precLatent = wAutoInvDiag + kron(trans(paramLatentAuto(j,0))*wparamLatentAuto,diagMatLev*trans(diagMatLev));

			// Calculate the variance-covariance matrix for the latentAuto variables
			mat varLatent = inv_sympd(precLatent);

			// Calculate the means for the new latentAuto variables
			mat meansLatent = varLatent*vectorise(diagMatLev*Yresid*wparamLatentAuto);

			// Calculate the new latentAuto variables
			latentAuto(j,0) = reshape(rmvnorm(1, meansLatent,varLatent),nAutoLev(j),nLatentAuto(j));
		}
	}

	return latentAuto;
}

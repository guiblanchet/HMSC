#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "euclDist.h"

using namespace arma; 
using namespace Rcpp;

// Calculates the Euclidean distance among all pairs of samples (rows) of a matrix.
// This function constructs a square matrix.
arma::mat euclDist(arma::mat& X) {
	
	// Count the number of rows of X
	int nRows = X.n_rows;
	
	// Transpose X for speed
	mat Xt = trans(X);
	
	// Result object
	mat res(nRows,nRows,fill::zeros);
	
	// loop over each row to calculate the Euclidean distance 
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nRows; j++) {
			if(i>=j){
				res(i,j) = sqrt(accu(square(Xt.col(i)-Xt.col(j))));
			}
		}
	}
	
	// Generate symmetric matrix from res by reflecting the lower triangle to the upper triangle
	res = symmatl(res);

	// return result object
	return res;
}

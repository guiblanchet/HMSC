#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rwishart.h"

using namespace arma;
using namespace Rcpp;

//' @title Random Sample from a Wishart Distribution
//' @description Creates a random wishart distribution when given degrees of freedom and a scale matrix. 
//' @param df An \code{int}, which gives the degrees of freedom of the Wishart.  (> 0)
//' @param S A \code{matrix} with dimensions m x m that provides Sigma, the covariance matrix. 
//' @return A \code{matrix} that is a Wishart distribution, aka the sample covariance matrix of a Multivariate Normal Distribution
//' @seealso \code{\link{riwishart}} 
//' @author James J Balamuta
//' @examples 
//' #Call with the following data:
//' rwishart(3, diag(2))
//' 
//' # Validation
//' set.seed(1337)
//' S = toeplitz((10:1)/10)
//' n = 10000
//' o = array(dim = c(10,10,n))
//' for(i in 1:n){
//' o[,,i] = rwishart(20, S)
//' }
//' mR = apply(o, 1:2, mean)
//' Va = 20*(S^2 + tcrossprod(diag(S)))
//' vR = apply(o, 1:2, var)
//' stopifnot(all.equal(vR, Va, tolerance = 1/16))
//' 
arma::mat rwishart(unsigned int df, const arma::mat& S){
  // Dimension of returned wishart
  unsigned int m = S.n_rows;
  
  // Z composition:
  // sqrt chisqs on diagonal
  // random normals below diagonal
  // misc above diagonal
  arma::mat Z(m,m);
  
  // Fill the diagonal
  for(unsigned int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }
  
  // Fill the lower matrix with random guesses
  for(unsigned int j = 0; j < m; j++){  
    for(unsigned int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1);
    }
  }
  
  // Lower triangle * chol decomp
  arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
  
  // Return random wishart
  return C.t()*C;
}
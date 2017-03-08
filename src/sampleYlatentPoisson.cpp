#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentPoisson.h"

using namespace arma ;
using namespace Rcpp ;
//' @title Sample the model response matrix
//'
//' @description Sample the model response matrix after the Poisson link function was applied. This function is meant to be used internally.
//'
//' @param Y The species community matrix to be modelled.
//' @param Ylatent Model site by species community matrix after the link function is applied.
//' @param EstModel Estimated model for the site by species community matrix.
//' @param residVar Vector of parameters with as many entries as there are species. Each values in this vector describes how much variation is not explained by the model. This vector should only contains 1s, but it was included in this function to deal with potential situations that may arise (in other words, in case I forgot a special case).
//' @param nsp Numeric. Number of species in \code{Y}.
//' @param nsite Numeric. Number of site sampled.
//'
//' @export
// [[Rcpp::export]]
arma::mat sampleYlatentPoisson(arma::mat& Y,
							  arma::mat& Ylatent,
							  arma::mat& EstModel,
							  arma::vec residVar,
							  double nsp,
							  int nsite){

	vec dsi2 = residVar;
	mat si2 = repmat(dsi2,1,nsite).t();

//	si2.elem(si2>100) = 100;
//	si2.elem(si2<0.001) = 0.001;

	mat Ez1 = EstModel;
	mat Y1 = Y;

	double r = 1000;
	double logr = log(r);

	//posterior w
	mat z = Ylatent;
	mat Psi = z - logr;
	int n1 = nsite;
	int n2 = nsp;
	double b = r;
	mat c = abs(Psi);
	mat m = b /2 /c % tanh(c/2);

	//when c is too big, this converges to 2;
	mat v_part2 = (sinh(c)-c) / square(cosh(c/2));

	uvec c700 = find(c > 700);
	vec two(sum(c700));
	two.fill(2);
	v_part2.elem(c700) = two; // approximate for when the above expression starts fail numerically to NaN
	mat v = b /4 /pow(c,3) % v_part2;
	mat w = randn(n1,n2) % sqrt(v) + m;
	w = abs(w);

	//posterior Psi
	mat Psi_var = 1 / (w + 1/si2);
	mat Psi_mean = Psi_var % ( Y1 - r/2.0 + (Ez1 - logr) / si2);
//	Psi = normrnd( Psi_mean, sqrt(Psi_var));
	Psi = randn(n1,n2) % sqrt(Psi_var) + Psi_mean;
	z =  Psi + logr;

	// Return Ylatent
	return z;
}

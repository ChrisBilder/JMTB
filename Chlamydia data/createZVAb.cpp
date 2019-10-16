#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

// See McMahan's MHS_GTR_IP file for supplemental information
// N = number of individuals
// np = number of pools current individual is in
// p = current vector of probabilities, i.e. inverse link function
// Y = matrix with pool information for individuals
// Z = matrix with observations and other information

NumericVector createZVAb(int N,int q,NumericVector ZVAb,NumericMatrix ZVA,NumericMatrix Y,NumericMatrix b) {	

	int id;

	for(int i=0; i<N; i++){
		id = Y(i,2) - 1;	// site id
		for(int j=0; j<q; j++){
			ZVAb(i) += (ZVA(i,j) * b(j,id));
		}
	}
	return ZVAb;
}
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

NumericMatrix createU(int N,int q,NumericMatrix Z,NumericMatrix Y,NumericMatrix V,NumericMatrix b,NumericMatrix U) {	

	int id, index;

	for(int i=0; i<N; i++){
		id = Y(i,2) - 1;	// site id
		index = 0;
		for(int j=0; j<(q-1); j++){
			for(int m=(j+1); m<q; m++){
				U(i,index) = b(j,id) * Z(i,m) * V(m,m);
				index += 1;
			}
		}
	}
	return U;
}
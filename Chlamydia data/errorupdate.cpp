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
// na = number of assays
// J = number of groups
// Y = matrix with pool information for individuals
// G = matrix with observations and other information

List errorupdates(int N,int K,NumericMatrix Y,NumericMatrix G,NumericMatrix PSe,NumericMatrix PSp,int na) {

	List res;
	int Ysum, Gk, ck, as, tmp;
	
	for(int k=0; k<K; k++){
		Ysum = 0;
		Gk = G(k,0);
		ck = G(k,1);
		as = G(k,2) - 1;	// assay used
		for(int i=0; i<ck; i++){
			tmp = G(k,i+3) - 1;	// indices of individuals
			Ysum += Y(tmp,0);
		}
		if(Ysum > 0){
			if(Gk > 0)
				PSe(as,0) += 1;	// for alpha term
			else
				PSe(as,1) += 1;	// for beta term
		}
		else{
			if(Gk > 0)
				PSp(as,1) += 1;	// for alpha term
			else
				PSp(as,0) += 1;	// for beta term
		}
	}

	res["PSe"] = PSe;
	res["PSp"] = PSp;
	return res;
}
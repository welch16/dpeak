// Calculate score function to obtain updated mu (Rcpp)

#include <Rcpp.h> 

RcppExport SEXP cpp_fraglen( SEXP Lgrid, SEXP S, SEXP E, SEXP strand, SEXP Z, SEXP mu, SEXP nGroup, SEXP R, SEXP gamma ) { 
	using namespace Rcpp;
	
   Rcpp::NumericVector LgridVec(Lgrid), SVec(S), EVec(E), strandVec(strand), muVec(mu);
   Rcpp::NumericMatrix zMat(Z);
   double nGroupValue = Rcpp::as<double>(nGroup);
   double rValue = Rcpp::as<double>(R);
   double gammaValue = Rcpp::as<double>(gamma);
   
   Rcpp::NumericVector out( LgridVec.size() );
   
   for (int i = 0; i < LgridVec.size(); i++) {
     out[i] = 0;
     for (int j = 0; j < SVec.size(); j++) {
       for (int g = 0; g < nGroupValue; g++) {
	       if ( strandVec[j] == 1 ) {
		       // forward strand
		       
		       if ( SVec[j] <= muVec[g] && muVec[g] <= SVec[j] + LgridVec[i] - 1 ) {
			       out[i] += zMat(j,g) * log( ( 1 - gammaValue ) / LgridVec[i] );
		       } else {
			       out[i] += zMat(j,g) * log( gammaValue / rValue );
		       }
	       } else {
		       // backward strand
		       
		       if ( EVec[j] - LgridVec[i] + 1 <= muVec[g] && muVec[g] <= EVec[j] ) {
			       out[i] += zMat(j,g) * log( ( 1 - gammaValue ) / LgridVec[i] );
		       } else {
			       out[i] += zMat(j,g) * log( gammaValue / rValue );
		       }
	       } // loop: if, strand
       } // loop: for, group
     } // loop: for, N
   } // loop: for, L grid
   
   return out;
}

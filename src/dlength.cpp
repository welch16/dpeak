// Calculate P(L) (Rcpp)

#include <Rcpp.h> 

RcppExport SEXP cpp_dlength( SEXP L, SEXP Lname, SEXP Lfreq ) { 
	using namespace Rcpp;
	
   Rcpp::NumericVector LVec(L), LnameVec(Lname), LfreqVec(Lfreq);
   double denom = 0;   
   Rcpp::NumericVector out( LVec.size() );
   
   // calculate denominator
   
   for (int i = 0; i < LfreqVec.size(); i++ ) {
	   denom += LfreqVec[i];
   }
   
   // calculate P(L)
   
   for (int i = 0; i < LVec.size(); i++) {     
     for (int j = 0; j < LnameVec.size(); j++) {
	     if ( LVec[i] == LnameVec[j] ) {
		     out[i] = LfreqVec[j] / denom;
		     break;
	     }
     }
   }
   
   return out;
}

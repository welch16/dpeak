// Calculate score function to obtain updated mu (Rcpp)

#include <Rcpp.h> 

RcppExport SEXP cpp_score( SEXP grid, SEXP S, SEXP E, SEXP L, SEXP Zg, SEXP R, SEXP gamma ) { 
	using namespace Rcpp;
	
   Rcpp::NumericVector gridVec(grid), SVec(S), EVec(E), LVec(L), ZgVec(Zg);
   Rcpp::NumericVector out( gridVec.size() );
   double rValue = Rcpp::as<double>(R);
   double gammaValue = Rcpp::as<double>(gamma);
   
   for (int i = 0; i < gridVec.size(); i++) {
     out[i] = 0;
     for (int j = 0; j < SVec.size(); j++) {
       if ( SVec[j] <= gridVec[i] && gridVec[i] <= EVec[j] ) {
         out[i] = out[i] + ZgVec[j] * log( ( 1 - gammaValue ) / LVec[j] );
       } else {
         out[i] = out[i] + ZgVec[j] * log( gammaValue / rValue );
       }
     }
   }
   return out;
}

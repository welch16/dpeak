// Is Z0 larger than each Z? (Rcpp)

#include <Rcpp.h> 

RcppExport SEXP cpp_ismaxZ0( SEXP Z0, SEXP Z ) { 
	using namespace Rcpp;
	
   Rcpp::NumericVector z0Vec( Z0 );
   Rcpp::NumericMatrix zMat( Z );   
   int rowComp = 0;   
   int nMax = 0;
   
   for (int i = 0; i < z0Vec.size(); i++) {     
     rowComp = 0;
     
     for (int j = 0; j < zMat.ncol(); j++) {
	     if ( z0Vec[i] > zMat(i,j) ) {
		     rowComp++;
	     }
     }
     
     if ( rowComp == zMat.ncol() ) {
	     nMax++;
     }
   }
   
   return Rcpp::wrap(nMax);
}

// Normalize each row (Rcpp)

#include <Rcpp.h> 

RcppExport SEXP cpp_normalize( SEXP Z ) { 
	using namespace Rcpp;
	
	Rcpp::NumericMatrix zMat( Z );   
	int N = zMat.nrow();
	int nGroup = zMat.ncol();	
	double rowsum;	
	Rcpp::NumericMatrix out( zMat.nrow(), zMat.ncol() );
	
	for (int i = 0; i < N; i++) {
		// row sum
		
		rowsum = 0;
		for (int g=0; g < nGroup; g++) {
			rowsum += zMat(i,g);
		}
		
		// normalize & random sampling
		
		if ( rowsum > 0 ) {
			for (int g=0; g < nGroup; g++) {
				out(i,g) = zMat(i,g) / rowsum;
			}
		} else {
			// if this row is all zero, return 0
			
			for (int g=0; g < nGroup; g++) {
				out(i,g) = 0;
			}
		}	
		
	}
	
	return out;
}

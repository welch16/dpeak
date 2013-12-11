// Calculate score function to obtain updated mu (Rcpp)

#include <Rcpp.h> 

RcppExport SEXP cpp_samp( SEXP Z ) { 
	using namespace Rcpp;
	
	Rcpp::NumericMatrix zMat( Z );   
	int N = zMat.nrow();
	int nGroup = zMat.ncol();	
	Rcpp::NumericVector cumsum( nGroup );	
	double rVal;
	Rcpp::NumericVector out( N );
	
	for (int i = 0; i < N; i++) {
		out[i] = 0;
		
		// cumulative sum
		
		cumsum[0] = zMat(i,0);
		for (int g=1; g < nGroup; g++) {
			cumsum[g] = cumsum[(g-1)] + zMat(i,g);
		}
		
		// normalize & random sampling
		
		if ( cumsum[nGroup-1] > 0 ) {
			for (int g=0; g < nGroup; g++) {
				cumsum[g] /= cumsum[nGroup-1];
			}	
			
			// random sampling
			
			rVal = ::unif_rand();
			for (int g=0; g < nGroup; g++) {
				if ( cumsum[g] < rVal ) {
					out[i]++;
				}
			}
		}
		
		// if this row is all zero, return 0
	}
	
	return out;
}

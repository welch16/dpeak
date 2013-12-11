# choose optimal model

.selectModel <- function( BIC_vec, maxComp, pConst=0.2 ) { 

    # remove models with NA
    
    gvec <- c(1:maxComp)[ !is.na(BIC_vec) ]
    glen <- length(gvec)
    BIC_vec <- BIC_vec[ !is.na(BIC_vec) ]
    
    # find the first local min
    
    if ( glen <= 2 ) {
        # if length is shorter than two, just use min
        
        min_model <- gvec[ which.min(BIC_vec) ]
    } else {        
        min_model <- gvec[1]
        
        diff_bic <- diff(BIC_vec)                
        for ( k in 2:(glen-1) ) {
            if ( diff_bic[k-1] < 0 & diff_bic[k] > 0 ) {
                min_model <- gvec[k]
                break
            }
        }
        
        if ( min_model == gvec[1] & which.min(BIC_vec) == glen ) {
            min_model <- gvec[glen]
        }
    }    
    
    # find the smallest model which is similar to the first local min
    
    if ( glen <= 2 | which( gvec == min_model ) == 1 ) {
        opt_model <- min_model
    } else {
        # cutoff to determine plateau
        
        cutoff <- ( max(BIC_vec) - min(BIC_vec) ) * pConst
        
        # search plateau
        
        loc_min <- which( gvec == min_model )
        
        k <- loc_min - 1
        while( k >= 1 ) {
            if ( abs( BIC_vec[loc_min] - BIC_vec[k] ) <= cutoff ) {            
                loc_min <- k
                k <- loc_min - 1
            }
            k <- k - 1
        }
        
        opt_model <- gvec[loc_min]
    }
    
    return( opt_model )
}

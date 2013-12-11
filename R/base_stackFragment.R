
# stack fragment (projection to coordinates)

.stackFragment <- function( fragmentMat ) {      
    
    if ( !is.na(fragmentMat[1,1]) ) {
        xvarq <- c( min(fragmentMat[,1]), max(fragmentMat[,2]) )
        xvar <- c(xvarq[1]:xvarq[2])
        yvar <- .ff_stack( fragmentMat[,1], fragmentMat[,2], xvarq[1], xvarq[2] )   
        
        outMat <- cbind( xvar, yvar )
    } else {
        outMat <- matrix( NA )
    }
    
    return( outMat )
}

# GOF plot

.plotBIC <- function( object ) {

    
    maxComp <- object@maxComp
    
    # generate matrix of BIC & AIC
    
    bicMat <- aicMat <- matrix( NA, length(object@fits), maxComp )
    
    for ( i in 1:length(object@fits) ) {
        if ( any(is.na(object@optMu[[i]])) ) {
            # error treatment: skip peaks with no fragments
        } else {        
            bicMat[i,] <- object@bicVec[[i]]
            aicMat[i,] <- object@aicVec[[i]]
        }
    }
    
    # GOF plot
    for ( i in 1:length(object@stackedFragment) ) {   
        # error treatment: skip peaks with no fragments
        
        # flag if there are no reads in the peak region
        
        if ( is.null(object@fragSet[[i]]) ) {        
            par( mfrow=c(1,1) )
            
            plot_title <- paste(object@peakChr[i],": ",
                object@peakStart[i],"-",object@peakEnd[i],sep="")        
            
            plot( 0, 0, type="n", xlab="", ylab="", axes=FALSE,
                main=plot_title, xlim=c(-5,5), ylim=c(-5,5) )    
            text( 0, 0, "No fragment (read) in the peak region" )
            
            next;
        }
        
        # plot
        
        par( mfrow=c(1,2) )
        
        plot_title <- paste("BIC\n",
            object@peakChr[i],": ",object@peakStart[i],"-",object@peakEnd[i],sep="")        
        plot( 1:maxComp, bicMat[i,],  
            xlab="Number of components", ylab="BIC value", main=plot_title )
        xsub <- c(1:maxComp)[ !is.na(bicMat[i,]) ]
        ysub <- bicMat[i,][ !is.na(bicMat[i,]) ]
        lines( xsub, ysub, type="l", lty=2 )
        
        plot_title <- paste("AIC\n",
            object@peakChr[i],": ",object@peakStart[i],"-",object@peakEnd[i],sep="")        
        plot( 1:maxComp, aicMat[i,], 
            xlab="Number of components", ylab="AIC value", main=plot_title )
        xsub <- c(1:maxComp)[ !is.na(aicMat[i,]) ]
        ysub <- aicMat[i,][ !is.na(aicMat[i,]) ]
        lines( xsub, ysub, type="l", lty=2 )
    }
}

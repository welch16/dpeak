# plot fitting results

.plotFit <- function( object, threshold=1000, 
    strand=FALSE, extension=1, smoothing=FALSE ) {       
    
    # extract estimates
    
    PET <- object@PET 
    psize <- object@psize
    #delta <- floor( psize / 2 )
        
    # error treatment
    
    if ( extension < 1 ) {
        stop( "Negative 'extension' is not allowed!" )
    }
    
    # plot
    
    for ( i in 1:length(object@stackedFragment) ) {    
    
        plot_title <- paste(object@peakChr[i],": ",object@peakStart[i],"-",object@peakEnd[i],sep="")
    
        # flag if there are no reads in the peak region
        
        if ( is.na(object@fragSet[[i]][1,1]) ) {
                
            plot( 0, 0, type="n", xlab="", ylab="", axes=FALSE,
                main=plot_title, xlim=c(-5,5), ylim=c(-5,5) )    
            text( 0, 0, "No fragment (read) in the peak region" )
            
            next;
        }
        
        # plot
            
        xlim <- rep( NA, 2 )
        xlim[1] <- min( object@peakStart[i], object@stackedFragment[[i]][,1] )
        xlim[2] <- max( object@peakEnd[i], object@stackedFragment[[i]][,1] )
        
        Ni <- nrow(object@fragSet[[i]])
            
        if ( strand==TRUE ) {
            
            .plotStrandData( stackedFragment=object@stackedFragment[[i]],
                fragSet=object@fragSet[[i]], plot_title=plot_title, xlim=xlim,
                PET=PET, extension=extension, smoothing=smoothing )
                
            if ( extension == 1 ) {
                legend( "topright", lty=c(1,1,2,2), lwd=c(1,1,2,2),
                    col=c("black","red","blue","lightblue"),
                    c("Stacked reads (forward)","Stacked reads (reverse)",
                        paste("Binding sites (strength > ",threshold,")",sep=""),
                        paste("Binding sites (strength <= ",threshold,")",sep=""))
                )
            } else {
                legend( "topright", lty=c(1,1,2,2), lwd=c(1,1,2,2),
                    col=c("black","red","blue","lightblue"),
                    c("Stacked extended reads (forward)","Stacked extended reads (reverse)",
                        paste("Binding sites (strength > ",threshold,")",sep=""),
                        paste("Binding sites (strength <= ",threshold,")",sep=""))
                )
            }
        } else {
            # combined
            
            plot( object@stackedFragment[[i]][,1], object@stackedFragment[[i]][,2], type="l", 
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title, 
                xlim=xlim, ylim=c(0,max(object@stackedFragment[[i]][,2])*1.2) )            
                
            legend( "topright", lty=c(1,2,2), lwd=c(1,2,2),
                col=c("black","blue","lightblue"), 
                c("Stacked fragments",
                    paste("Binding sites (strength > ",threshold,")",sep=""),
                    paste("Binding sites (strength <= ",threshold,")",sep="")
                )
            )
        }            

        abline( v=object@optMu[[i]][ Ni * object@optPi[[i]] > threshold ], col="blue", lty=2, lwd=2 )
        abline( v=object@optMu[[i]][ Ni * object@optPi[[i]] <= threshold ], col="lightblue", lty=2, lwd=2 )
        
        #for ( j in 1:length(object@optMu[[i]]) ) {
        #    if ( object@optPi[[i]][j] > threshold ) {
        #        rect( object@optMu[[i]][j] - delta, -10, 
        #            object@optMu[[i]][j] + delta, max(object@stackedFragment[[i]][,2])+20, 
        #            density=10, col="blue", border="blue" )
        #    } else {
        #        rect( object@optMu[[i]][j] - delta, -10, 
        #            object@optMu[[i]][j] + delta, max(object@stackedFragment[[i]][,2])+20, 
        #            density=10, col="lightblue", border="lightblue" )
        #    }
        #}
        
        abline( v=object@peakStart[i], col="red", lty=2 )
        abline( v=object@peakEnd[i], col="red", lty=2 )  
    }
}

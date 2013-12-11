# GOF plot

.plotGOF <- function( object, nsimul=10000, seed=12345, 
    extension=1, smoothing=TRUE, nCore=8 ) {
    
    # extract estimates
    
    maxComp <- object@maxComp
    #bgComp <- object@bgComp
    peakStart <- object@peakStart
    peakEnd <- object@peakEnd
    
    PET <- object@PET
    fragLenTable <- object@fragLenTable
    aveFragLen <- object@aveFragLen
    Fratio <- object@Fratio
    
    # generate list to generate fragments
    
    optList <- vector( "list", length(object@optMu) )
    for ( i in 1:length(object@optMu) ) {
        # error treatment: skip peaks with no fragments
        
        if ( is.na(object@fragSet[[i]][1,1]) ) {
            optList[[i]] <- matrix( NA )
        } else {        
            # stack info
            
            optList_i <- list()
            
            optList_i$nsimul <- nsimul
            optList_i$mu <- object@optMu[[i]]
            optList_i$pi <- object@optPi[[i]]
            optList_i$pi0 <- object@optPi0[[i]]
            optList_i$minS <- min(object@stackedFragment[[i]][,1])
            optList_i$maxS <- max(object@stackedFragment[[i]][,1])
            optList_i$peakstart <- peakStart[i]
            optList_i$peakend <- peakEnd[i]
            if ( PET == FALSE ) {
                optList_i$delta <- object@optDelta[[i]]
                optList_i$sigma <- object@optSigma[[i]]
            }            
            optList[[i]] <- optList_i
        }
    }
    
    if ( PET ) {
        Lvalue <- as.numeric(as.vector(names(fragLenTable)))
        Lprob <- as.numeric(fragLenTable)
    }
    
    # generate fragments (using parallel computing, if parallel exists)
    
    message( "Info: Generating simulated fragments from the fitted model..." )
    
    if ( is.element( "parallel", installed.packages()[,1] ) ) {
        # if "parallel" package exists, utilize parallel computing with "mclapply"
        library(parallel)
        
        simList <- mclapply( optList, function(x) { 
            simFrag <- .generateFragment( object=x, seed=seed, 
                PET=PET, Lvalue=Lvalue, Lprob=Lprob, 
                Fratio=Fratio, aveFragLen=aveFragLen ) 
            simStack <- .stackFragment( simFrag )
            out <- list( stackedFragment=simStack, fragSet=simFrag )
            return( out )
        }, mc.cores = nCore )
    } else {
        # otherwise, use usual "lapply"
        
        simList <- lapply( optList, function(x) { 
            #x <- optList[[i]]
            simFrag <- .generateFragment( object=x, seed=seed, 
                PET=PET, Lvalue=Lvalue, Lprob=Lprob, 
                Fratio=Fratio, aveFragLen=aveFragLen ) 
            simStack <- .stackFragment( simFrag )
            out <- list( stackedFragment=simStack, fragSet=simFrag )
            return( out )
        } )
    }    
    
    # GOF plot
    
    message( "Info: Generating GOF plots..." )
    
    for ( i in 1:length(object@stackedFragment) ) {   
            
        plot_title <- paste(object@peakChr[i],": ",
            object@peakStart[i],"-",object@peakEnd[i],sep="")
        
        # flag if there are no reads in the peak region
        
        if ( is.na(object@fragSet[[i]][1,1]) ) {
                
            plot( 0, 0, type="n", xlab="", ylab="", axes=FALSE,
                main=plot_title, xlim=c(-5,5), ylim=c(-5,5) )    
            text( 0, 0, "No fragment (read) in the peak region" )
            
            next;
        }
            
        # plot
        
        stackedSimFrag <- simList[[i]]$stackedFragment
        
        xlim <- rep( NA, 2 )
        xlim[1] <- min( object@peakStart[i], object@stackedFragment[[i]][,1], 
            stackedSimFrag[,1] )
        xlim[2] <- max( object@peakEnd[i], object@stackedFragment[[i]][,1], 
            stackedSimFrag[,1] )
        
        if ( PET == FALSE ) {        
            .plotStrandData( stackedFragment=object@stackedFragment[[i]],
                fragSet=object@fragSet[[i]], plot_title=plot_title, xlim=xlim,
                PET=PET, extension=extension, smoothing=smoothing )
                
            #if ( extension == 1 ) {
            #    legend( "topright", lty=c(1,1), col=c("black","red"),
            #        c("Stacked reads (forward)","Stacked reads (reverse)")
            #    )
            #} else {
            #    legend( "topright", lty=c(1,1), col=c("black","red"),
            #        c("Stacked extended reads (forward)","Stacked extended reads (reverse)")
            #    )                
            #}
        } else {
            # PET
            
            plot( object@stackedFragment[[i]][,1], 
                object@stackedFragment[[i]][,2], type="l", 
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title, 
                xlim=xlim, ylim=c(0,max(object@stackedFragment[[i]][,2])*1.2) )            
                
            #legend( "topright", lty=1, col="black", "Stacked fragments" )
        }            
        abline( v=object@peakStart[i], col="red", lty=2 )
        abline( v=object@peakEnd[i], col="red", lty=2 )  
        
        # plot simulated data
                
        if ( PET == TRUE ) {
            # PET
            
            simX <- stackedSimFrag[,1]
            simY <- stackedSimFrag[,2] * max(object@stackedFragment[[i]][,2]) /
                max(stackedSimFrag[,2])
                # adjust frequency to make frequency comparable to observed one
            lines( simX, simY, col="gray", lty=2 )
            
            #abline( v=peakStart[i], col="red", lty=2 )
            #abline( v=peakEnd[i], col="red", lty=2 )  
                
            legend( "topright", lty=c(1,2), col=c("black","gray"),
                c("Stacked fragments (data)","Stacked fragments (simulated)") )
        } else {
            # SET
            
            ratio <- max(object@stackedFragment[[i]][,2]) /
                max(stackedSimFrag[,2])
                  
            .lineStrandData( stackedFragment=stackedSimFrag,
                fragSet=simList[[i]]$fragSet, adjustment=ratio,
                plot_title=plot_title, xlim=xlim,
                PET=PET, extension=extension, smoothing=smoothing )
                
            if ( extension == 1 ) {
                legend( "topright", lty=c(1,1,2,2), 
                    col=c("black","red","gray","pink"),
                    c("Stacked reads: forward (data)",
                    "Stacked reads: reverse (data)",
                    "Stacked reads: forward (simulated)",
                    "Stacked reads: reverse (simulated)")
                )
            } else {
                legend( "topright", lty=c(1,1,2,2), 
                    col=c("black","red","gray","pink"),
                    c("Stacked extended reads: forward (data)",
                    "Stacked extended reads: reverse (data)",
                    "Stacked extended reads: forward (simulated)",
                    "Stacked extended reads: reverse (simulated)")
                )
            }
        }
    }
    
    message( "Info: Done!" )
}

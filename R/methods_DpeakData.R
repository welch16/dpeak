
# generic methods for "DpeakData" class

setMethod(
    f="show",
    signature="DpeakData",
    definition=function( object ) {    
        # extract objects
        
        peakChr <- object@peakChr
        chrList <- unique(peakChr)
        #fragSet <- object@fragSet
        PET <- object@PET
        aveFragLen <- object@aveFragLen
        
        # summary
        nFrag <- unlist( lapply( object@fragSet, 
            function(x) ifelse( !is.null(x), length(x), 0 )
         ) )
        sumRead <- sum( nFrag )
        medNumRead <- median( nFrag )
        
        cat( "------------------------------------------------------------\n" )
        cat( "Summary: Dpeak data (class: DpeakData)\n" )
        cat( "------------------------------------------------------------\n" )
        cat( "Number of peaks: ",length(peakChr),"\n", sep="" )
        cat( "Number of chromosomes: ",length(chrList),"\n", sep="" )
        cat( "Tag type: ",ifelse(PET,"PET","SET"),"\n", sep="" )
        if ( PET == TRUE ) {
            cat( "Median fragment length: ",aveFragLen,"\n", sep="" )    
        } else {
            cat( "Fragment length (provided by user): ",aveFragLen,"\n", sep="" )    
        }
        cat( "Number of utilized reads: ",sumRead,"\n", sep="" )
        cat( "Median number of reads in each peak: ",medNumRead,"\n", sep="" )
        cat( "------------------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="DpeakData",
    definition=function( x ) {
        warning( "'print' method for 'DpeakData' class is not supported yet." )
    }
)

setMethod(
    f="printEmpty",
    signature="DpeakData",
    definition=function( object, ... ) {
        if ( object@emptyList[1] == "" ) {
            message( "Every peak region has at least one read." )
        } else {
            out <- apply( as.matrix(object@emptyList), 1, 
                function(x) { 
                    splitOrg <- strsplit( x, "_" )[[1]]
                    
                    if ( length(splitOrg) > 3 ) {
                        xSplit <- rep( NA, 3 )
                        xSplit[1] <- paste( splitOrg[ 1:(length(splitOrg)-2) ],
                            collapse="_" )
                        xSplit[2:3] <- splitOrg[ -c(1:(length(splitOrg)-2)) ]
                    } else {
                        xSplit <- splitOrg
                    }
                    return(xSplit)
            } )
            
            out <- data.frame( t(out), stringsAsFactors=FALSE )
            out[,2] <- as.numeric(out[,2])
            out[,3] <- as.numeric(out[,3])
            colnames(out) <- c( "chrID", "peakStart", "peakEnd" )
            
            return( out )
        }
    }
)

setMethod(
    f="plot",
    signature=c("DpeakData","missing"),
    definition=function( x, y, filename=NULL, strand=FALSE, extension=1, smoothing=FALSE, ... ) {        
        # extract estimates
        
        PET <- x@PET    
        
        # error treatment
        
        if ( extension < 1 ) {
            stop( "Negative 'extension' is not allowed!" )
        }
        
        # plot
        
        pdf( filename ) 
        
        for ( i in 1:length(x@stackedFragment) ) {   
        
            plot_title <- paste(x@peakChr[i],": ",x@peakStart[i],"-",x@peakEnd[i],sep="")
            
            # flag if there are no reads in the peak region
            
            if ( is.na(x@fragSet[[i]][1,1]) ) {
                    
                plot( 0, 0, type="n", xlab="", ylab="", axes=FALSE,
                    main=plot_title, xlim=c(-5,5), ylim=c(-5,5) )    
                text( 0, 0, "No fragment (read) in the peak region" )
                
                next;
            }
            
            # plot
            
            xlim <- rep( NA, 2 )
            xlim[1] <- min( x@peakStart[i], x@stackedFragment[[i]][,1] )
            xlim[2] <- max( x@peakEnd[i], x@stackedFragment[[i]][,1] )
            
            if ( strand==TRUE ) {
            
                .plotStrandData( stackedFragment=x@stackedFragment[[i]],
                    fragSet=x@fragSet[[i]], plot_title=plot_title, xlim=xlim,
                    PET=PET, extension=extension, smoothing=smoothing )
                    
                if ( extension == 1 ) {
                    legend( "topright", lty=c(1,1), col=c("black","red"),
                        c("Stacked reads (forward)","Stacked reads (reverse)")
                    )
                } else {
                    legend( "topright", lty=c(1,1), col=c("black","red"),
                        c("Stacked extended reads (forward)","Stacked extended reads (reverse)")
                    )                
                }
            } else {
                # combined
                
                plot( x@stackedFragment[[i]][,1], x@stackedFragment[[i]][,2], type="l", 
                    xlab="Genomic coordinates", ylab="Frequency",
                    main=plot_title, 
                    xlim=xlim, ylim=c(0,max(x@stackedFragment[[i]][,2])*1.2) )            
                    
                legend( "topright", lty=1, col="black", "Stacked fragments" )
            }            
            abline( v=x@peakStart[i], col="red", lty=2 )
            abline( v=x@peakEnd[i], col="red", lty=2 )  
        }
        
        dev.off()        
    }
)

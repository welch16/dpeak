# plot fitting results

.plotStrandData <- function( stackedFragment, fragSet, plot_title, xlim,
    PET, extension=1, smoothing=FALSE ) {

   xvar <- xlim[1]:xlim[2]
   ## xvar <- stackedFragment[,1]
    

   ## xval <- cumsum(stackedFragment[[3]][[1]])  + stackedFragment[[1]] -1
   ##          yval <- stackedFragment[[3]][[2]]
   ##          if(xlim[1] < min(xval)){
   ##            xval <- c(xlim[1],xval)
   ##            yval <- c(0,yval)
   ##          }
   ##          if(xlim[2] > max(xval)){
   ##            xval <- c(xval,xlim[2])
   ##            yval <- c(yval,0)
   ##          }
   ##          x <- xlim[1]:xlim[2]
   ##          y <- stepfun(xval,c(yval,0))(x)

    
    # forward strand

    if ( PET == FALSE ) {
        # SET
        
        ## indF <- which( fragSet[,3] == "F" )
        indF <- which(strand(fragSet) == "+")
        
        if ( length(indF) > 0 ) {                
            ## Fstart <- fragSet[ indF, 1 ]
            ## Fend <- pmin( fragSet[ indF, 1 ] + extension - 1, max(xvar) )
          Fstart <- start(fragSet)[indF]
          Fend <- pmin(Fstart + extension -1,max(xvar))
                
            xvarq <- c( min(Fstart), max(Fend) )
            xvar_org <- c(xvarq[1]:xvarq[2])
            yvar_org <- .ff_stack( Fstart, Fend, xvarq[1], xvarq[2] ) 
            
            yvarF <- rep( 0, length(xvar) )
            yvarF[ match( xvar_org, xvar ) ] <- yvar_org
        } else {
            yvarF <- c()
        }
    } else {
        # PET
      #### need to fix this  
        Fstart <- start(fragSet)
        Fend <- pmin( start(fragSet) + extension - 1, max(xvar) )
                
        xvarq <- c( min(Fstart), max(Fend) )
        xvar_org <- c(xvarq[1]:xvarq[2])
        yvar_org <- .ff_stack( Fstart, Fend, xvarq[1], xvarq[2] ) 
        
        yvarF <- rep( 0, length(xvar) )
        yvarF[ match( xvar_org, xvar ) ] <- yvar_org
    }
    
    # reverse strand
        
    if ( PET == FALSE ) {
        # SET
        
        if ( length(fragSet) - length(indF) > 0 ) {
            if ( length(indF) > 0 ) {
                ## Rstart <- pmax( fragSet[ -indF, 2 ] - extension + 1, min(xvar) )
                ## Rend <- fragSet[ -indF, 2 ]
                Rstart <- pmax( end(fragSet)[ -indF ] - extension + 1, min(xvar) )
                Rend <- end(fragSet)[-indF]       
            } else {
                Rstart <- pmax( end(fragSet) - extension + 1, min(xvar) )
                Rend <- end(fragSet)
            }
                
            xvarq <- c( min(Rstart), max(Rend) )
            xvar_org <- c(xvarq[1]:xvarq[2])
            yvar_org <- .ff_stack( Rstart, Rend, xvarq[1], xvarq[2] ) 
            
            yvarR <- rep( 0, length(xvar) )
            yvarR[ match( xvar_org, xvar ) ] <- yvar_org
        } else {
            yvarR <- c()
        }
    } else {
        # PET
        #### this too
        Rstart <- pmax( end(fragSet)- extension + 1, min(xvar) )
        Rend <- end(fragSet)
                
        xvarq <- c( min(Rstart), max(Rend) )
        xvar_org <- c(xvarq[1]:xvarq[2])
        yvar_org <- .ff_stack( Rstart, Rend, xvarq[1], xvarq[2] ) 
        
        yvarR <- rep( 0, length(xvar) )
        yvarR[ match( xvar_org, xvar ) ] <- yvar_org
    }
    
    # plot
        
    if ( smoothing ) {
        if ( length(yvarF) > 0 ) {
            smoothF <- smooth.spline( xvar, yvarF )
        }
        if ( length(yvarR) > 0 ) {
            smoothR <- smooth.spline( xvar, yvarR )
        }
        
        if ( length(yvarF) > 0 & length(yvarR) > 0 ) {
            plot( smoothF, 
                type="l", col="black",
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title, xlim=xlim, ylim=c(0,max(smoothF$y,smoothR$y)*1.2) )
            lines( smoothR, col="red" )
        } else if ( length(yvarF) > 0 ) {
            plot( smoothF, 
                type="l", col="black",
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title, xlim=xlim, ylim=c(0,max(smoothF$y)*1.2) )
        } else if ( length(yvarR) > 0 ) {
            plot( smoothR, 
                type="l", col="red",
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title, xlim=xlim, ylim=c(0,max(smoothR$y)*1.2) )
        }
    } else {
        if ( length(yvarF) > 0 & length(yvarR) > 0 ) {
            plot( xvar, yvarF, 
                type="l", col="black",
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title, xlim=xlim, ylim=c(0,max(yvarF,yvarR)*1.2) )                    
            lines( xvar, yvarR, col="red" )
        } else if ( length(yvarF) > 0 ) {
            plot( xvar, yvarF, 
                type="l", col="black",
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title, xlim=xlim, ylim=c(0,max(yvarF)*1.2) )
        } else if ( length(yvarR) > 0 ) {
            plot( xvar, yvarR,
                type="l", col="red",
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title, xlim=xlim, ylim=c(0,max(yvarR)*1.2) )
        }
    }       
}


# generic methods for "DpeakFit" class

setMethod(
    f="show",
    signature="DpeakFit",
    definition=function( object ) {    
        # extract objects
        
        optMu <- object@optMu
        optPi0 <- object@optPi0
        maxComp <- object@maxComp
        #bgComp <- object@bgComp
        
        # summary
        
        med_opt_mu <- median( unlist( lapply( optMu, 
            function(x) ifelse( any(!is.na(x)), length(x), 0 )
        ) ) )
        #nonzero_pi0_list <- unlist( lapply( optPi0, 
        #    function(x) ifelse( any(!is.na(x)), as.numeric(x>0), 0 )
        #) )
        #nonzero_pi0 <- round( 100 * sum(nonzero_pi0_list) / length(nonzero_pi0_list) ) 
        pi0_vec <- unlist( optPi0 )
        med_explain <- round( 100 * median( 1 - pi0_vec, na.rm=TRUE ) ) / 100
        
        cat( "------------------------------------------------------------\n" )
        cat( "Summary: Dpeak model fitting (class: DpeakFit)\n" )
        cat( "------------------------------------------------------------\n" )
        cat( "- Maximum possible number of binding events in each peak: ",maxComp,"\n", sep="" )
        cat( "- Median number of binding events in each peak: ",med_opt_mu,"\n", sep="" )
        #cat( "- Use background components to estimate binding sites? ",ifelse(bgComp,"Yes","No"),"\n", sep="" )
        #cat( "- Percentage of peaks with non-zero background proportion: ",nonzero_pi0," %\n", sep="" )
        cat( "- Median explanation ratio: ",med_explain," %\n", sep="" )
        cat( "------------------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="DpeakFit",
    definition=function( x ) {
        warning( "'print' method for 'DpeakFit' class is not supported yet." )
    }
)

setMethod(
    f="plot",
    signature=c("DpeakFit","missing"),
    definition=function( x, y, filename=NULL, plotType="fit", 
        strand=FALSE, extension=1, smoothing=FALSE,
        threshold=1000, nsimul=10000, seed=12345, nCore=8, ... ) {      

      browser()
      
        pdf( filename )
        
        if ( plotType == "fit" ) {
            
            # fitting results
            
            .plotFit( object=x, threshold=threshold, 
                strand=strand, extension=extension, smoothing=smoothing )
        } else if ( plotType == "GOF" ) {  
        
            # error treatment
            
            #if ( plotType == "GOF" & strand == TRUE ) {
            #    stop( "'strand=TRUE' is not supported for GOF plot!" )
            #}
        
            # GOF plot
            
            .plotGOF( object=x, nsimul=nsimul, seed=seed, 
                extension=extension, smoothing=smoothing, nCore=nCore )
        } else if ( plotType == "BIC" ) {          
            # plot of BIC curve
            
            .plotBIC( object=x )
        } else {
            stop( "Inappropriate 'plotType'! Use either 'fit', 'GOF', or 'BIC'!" )
        }
        
        dev.off()        
    }
)

setMethod(
    f="export",
    signature="DpeakFit",
    definition=function( object, type=NA, filename=NA, ... ) {
        # error treatment: check invalid type
         
        if ( is.na(type) )
        {
            message( "Info: 'type' is not specified by the user." )
            message( "Info: 'type' is specified as 'bed' instead." )
            type <- "bed"        
        }
        
        allType <- c("txt","gff","bed")
        invalidType <- TRUE
        for ( i in 1:length(type) )
        {
            if ( length(which(!is.na(match(type[i],allType)))) == 0 )
            {
                invalidType <- FALSE
            }
        }
        if ( !invalidType )
        {
            message( "Info: 'type' incorrect." )
            message( "Info: 'type' is specified as 'bed' instead." )
            type <- "bed"
        }
        
        # error treatment: 'filename' not specified
        
        if ( is.na(filename) )
        {
            message( "Info: 'filename' is not specified by the user." )
            message( "Info: 'filename' is specified as 'bindingSiteList' instead." )
            filename <- paste("bindingSiteList.",type,sep="")
        }
        
        # extract objects
        
        optMu <- object@optMu
        optPi <- object@optPi
        peakChr <- object@peakChr
        peakStart <- object@peakStart
        peakEnd <- object@peakEnd
        psize <- object@psize
        
        # process fitting results for exporting
        
        peakName <- paste( peakChr, ":", peakStart, "-", peakEnd, sep="" )
        #chrVec <- muVec <- piVec <- nameVec <- c()
        chrVec <- muVec <- strengthVec <- nameVec <- c()
        for ( i in 1:length(optMu) ) {
            # error treatment: skip peaks with no fragments
            
            #if ( any(is.na(optMu[[i]])) ) {
            ## if ( is.na(object@fragSet[[i]][1,1]) ) {
            if(is.null(object@fragSet[[i]])){
                next;
            }
            
            # stack data
            
            chrVec <- c( chrVec, rep( peakChr[i], length(optMu[[i]]) ) )
            muVec <- c( muVec, optMu[[i]] )
            #piVec <- c( piVec, 1000 * optPi[[i]] )
            ## strengthVec <- c( strengthVec, nrow(object@fragSet[[i]]) * optPi[[i]] )
            strengthVec <- c( strengthVec, length(object@fragSet[[i]]) * optPi[[i]] )
            nameVec <- c( nameVec, rep( peakName[i], length(optMu[[i]]) ) )
        }

        # export peak lists 
        
        message( "Info: Exporting the binding site list..." )
        
        switch( type,
            "gff" = {
                #.exportGFF( chrVec=chrVec, muVec=muVec, piVec=piVec, psize=psize,
                .exportGFF( chrVec=chrVec, muVec=muVec, 
                    strengthVec=strengthVec, psize=psize,
                    nameVec=nameVec, filename=filename )
            },
            "bed" = {
                .exportBED( chrVec=chrVec, muVec=muVec, 
                    strengthVec=strengthVec, psize=psize,
                    nameVec=nameVec, filename=filename )
            },
            "txt" = {
                .exportTXT( chrVec=chrVec, muVec=muVec, 
                    strengthVec=strengthVec, psize=psize,
                    nameVec=nameVec, filename=filename )
            }
        )
    }
)

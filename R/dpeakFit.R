# user interface for deconvolution
# (support parallel computing using "parallel" package)

setMethod(
    f="dpeakFit",
    signature="DpeakData",
    definition=function( object, objectMotif=NULL, estDeltaSigma="common", init="localmax",
		nTop=100, lbDelta=25, lbSigma=25,
        psize=21, maxComp=5, pConst=0.2, 
        nCore=8, verbose=FALSE, seed=12345, iterInit=50, iterMain=25, epsilon=1e-6 )
    {    
		# use motif information for initialization?
		
		if ( !is.null(objectMotif) ) {
			message( "Info: positions of binding events are initialized using sequence information." )
		}
		
		# how to initialize binding events
		
		if ( init == "localmax" ) {
			message( "Info: positions of binding events are initialized based on signal strength." )
		} else if ( init == "uniform" ) {
			message( "Info: positions of binding events are initialized uniformly over the candidate region." )
		} else {
			stop( "Inappropriate 'init' argument! It should be either 'localmax' or 'uniform'." )
		}
		
        # safe guard: iterInit, iterMain
        
        if ( iterInit < 2 ) {
            message( "Info: 'iterInit' should be larger than or equal to 2. 'iterInit' is set to 2." )
            iterInit <- 2
        }
        if ( iterMain < 2 ) {
            message( "Info: 'iterMain' should be larger than or equal to 2. 'iterMain' is set to 2." )
            iterMain <- 2
        }
        
        # safe guard: estDeltaSigma
        
        if ( object@PET == FALSE ) {
            if ( estDeltaSigma == "separate" ) {
                message( "Info: estimate peak shape for each candidate region, separately." )
            } else if ( estDeltaSigma == "common" ) {
                message( "Info: estimate common peak shape using top candidate regions." )
            } else {
                stop( "Inappropriate 'estDeltaSigma' argument!" )
            }
        }
        
        # safe guard: lbDelta
        
        if ( length(lbDelta) == 1 ) {
	        lbDelta <- rep( lbDelta, 2 )
        } else if ( length(lbDelta) > 2 ) {
	        message( "Info: Length of 'lbDelta' can be at most 2. Only the first two numbers are used." )
	        lbDelta <- lbDelta[1:2]
        }
        
        if ( any( lbDelta < 5 ) ) {
            message( "Info: 'lbDelta' should be larger than or equal to 5. 'lbDelta' is set to 5." )
            lbDelta[ lbDelta < 5 ] <- 5
        }
        
        # safe guard: lbSigma
        
        if ( length(lbSigma) == 1 ) {
	        lbSigma <- rep( lbSigma, 2 )
        } else if ( length(lbSigma) > 2 ) {
	        message( "Info: Length of 'lbSigma' can be at most 2. Only the first two numbers are used." )
	        lbSigma <- lbSigma[1:2]
        }
        
        if ( any( lbSigma < 5 ) ) {
            message( "Info: 'lbSigma' should be larger than or equal to 5. 'lbSigma' is set to 5." )
            lbSigma[ lbSigma < 5 ] <- 5
        }
        
        # extract objects
        
        PET <- object@PET
        if ( PET == TRUE ) {
            L_table <- object@fragLenTable
            aveFragLen <- NA
        } else {
            L_table <- NA
            aveFragLen <- object@aveFragLen
        }
        Fratio <- object@Fratio
		
		# estimation of common peak shape
		
		if ( estDeltaSigma == "common" ) {
		
			# choose top candidate regions
			
			nread <- sapply( object@fragSet, nrow )
			nTopFinal <- min( nTop, length(nread) )
			nreadCutoff <- sort( nread, decreasing=TRUE )[ nTopFinal ]
        
			dataObj <- vector( "list", nTopFinal )
			selvec <- which( nread >= nreadCutoff )
			for ( i in 1:length(selvec) ) {
				isel <- selvec[i]
				
				dataObj[[i]] <- list()
				dataObj[[i]]$frag <- object@fragSet[[isel]]
				dataObj[[i]]$peak <- c( object@peakStart[isel], object@peakEnd[isel] )
				dataObj[[i]]$signal <- object@stackedFragment[[isel]]
			
				if ( !is.null(objectMotif) ) {
					dataObj[[i]]$locmotif <- objectMotif@locMotif[[isel]]
				} else {
					dataObj[[i]]$locmotif <- NA
				}
			}
			
			# deconvolve top candidate regions (using parallel computing, if parallel exists)
			
			if ( is.element( "parallel", installed.packages()[,1] ) ) {
				# if "parallel" package exists, utilize parallel computing with "mclapply"
				library(parallel)
				
				fit_top <- mclapply( dataObj, function(x) { 
					.deconWrapper( fData=x, estDeltaSigma="separate", init=init, 
						deltaInit=NA, sigmaInit=NA, lbDelta=lbDelta, lbSigma=lbSigma,
						psize=psize, max_comp=maxComp, pConst=pConst,
						seed=seed, niter_init=iterInit, niter_gen=iterMain, 
						PET=PET, L_table=L_table, Fratio=Fratio, aveFragLen=aveFragLen,
						stop_eps=epsilon, verbose=verbose ) 
					}, mc.cores = nCore )
			} else {
				# otherwise, use usual "lapply"
				
				fit_top <- lapply( dataObj, function(x) { 
					.deconWrapper( fData=x, estDeltaSigma="separate", init=init,  
						deltaInit=NA, sigmaInit=NA, lbDelta=lbDelta, lbSigma=lbSigma,
						psize=psize, max_comp=maxComp, pConst=pConst,
						seed=seed, niter_init=iterInit, niter_gen=iterMain, 
						PET=PET, L_table=L_table, Fratio=Fratio, aveFragLen=aveFragLen,
						stop_eps=epsilon, verbose=verbose ) 
					} )
			}    
			
			# estimate common delta & sigma
			
			deltaCommon <- sigmaCommon <- ntotal <- 0
			
			for ( itop in 1:nTopFinal ) {
				ni <- nrow(dataObj[[itop]]$frag)
				
				deltaCommon <- deltaCommon + fit_top[[itop]]$optDelta * ni
				sigmaCommon <- sigmaCommon + fit_top[[itop]]$optSigma * ni	
				ntotal <- ntotal + ni
			}
			
			deltaCommon <- deltaCommon / ntotal
			sigmaCommon <- sigmaCommon / ntotal
			
			rm( dataObj )
			gc()
			
        } else if ( estDeltaSigma == "separate" ) {
			
			deltaCommon <- NA
			sigmaCommon <- NA
		}
		
        # construct object for model fitting
        
        dataObj <- vector( "list", length(object@fragSet) )
        for ( i in 1:length(object@fragSet) ) {
            dataObj[[i]] <- list()
            dataObj[[i]]$frag <- object@fragSet[[i]]
            dataObj[[i]]$peak <- c( object@peakStart[i], object@peakEnd[i] )
			dataObj[[i]]$signal <- object@stackedFragment[[i]]
			
			if ( !is.null(objectMotif) ) {
				dataObj[[i]]$locmotif <- objectMotif@locMotif[[i]]
			} else {
				dataObj[[i]]$locmotif <- NA
			}
        }
        
        # deconvolve all peaks (using parallel computing, if parallel exists)
        
        if ( is.element( "parallel", installed.packages()[,1] ) ) {
            # if "parallel" package exists, utilize parallel computing with "mclapply"
            library(parallel)
            
            fit_all <- mclapply( dataObj, function(x) { 
                .deconWrapper( fData=x, estDeltaSigma=estDeltaSigma, init=init,  
					deltaInit=deltaCommon, sigmaInit=sigmaCommon, lbDelta=lbDelta, lbSigma=lbSigma,
                    psize=psize, max_comp=maxComp, pConst=pConst,
                    seed=seed, niter_init=iterInit, niter_gen=iterMain, 
                    PET=PET, L_table=L_table, Fratio=Fratio, aveFragLen=aveFragLen,
                    stop_eps=epsilon, verbose=verbose ) 
                }, mc.cores = nCore )
        } else {
            # otherwise, use usual "lapply"
            
            fit_all <- lapply( dataObj, function(x) { 
                .deconWrapper( fData=x, estDeltaSigma=estDeltaSigma, init=init,  
					deltaInit=deltaCommon, sigmaInit=sigmaCommon, lbDelta=lbDelta, lbSigma=lbSigma,
                    psize=psize, max_comp=maxComp, pConst=pConst,
                    seed=seed, niter_init=iterInit, niter_gen=iterMain, 
                    PET=PET, L_table=L_table, Fratio=Fratio, aveFragLen=aveFragLen,
                    stop_eps=epsilon, verbose=verbose ) 
                } )
        }    
        
        rm( dataObj )
        gc()
        
        # select optimal model
        
        fits <- vector( "list", length(object@fragSet) )
        optFit <- optMu <- optPi <- optPi0 <- 
            optGamma <- optDelta <- optSigma <- bicVec <- aicVec <-
            vector( "list", length(object@fragSet) )
        
        for ( i in 1:length(object@fragSet) ) {
            fits[[i]] <- fit_all[[i]]$fits
            optFit[[i]] <- fit_all[[i]]$optFit
            optMu[[i]] <- fit_all[[i]]$optMu
            optPi[[i]] <- fit_all[[i]]$optPi
            optPi0[[i]] <- fit_all[[i]]$optPi0
            if ( PET == FALSE ) {
                optDelta[[i]] <- fit_all[[i]]$optDelta
                optSigma[[i]] <- fit_all[[i]]$optSigma
            } else {
                optGamma[[i]] <- fit_all[[i]]$optGamma
            }
            bicVec[[i]] <- fit_all[[i]]$bicVec
            aicVec[[i]] <- fit_all[[i]]$aicVec
        }
        rm( fit_all )
        gc()
        
        names(fits) <- names(optFit) <- names(optMu) <- names(optPi) <- names(optPi0) <- 
        names(optGamma) <- names(optDelta) <- names(optSigma) <- names(bicVec) <- names(aicVec) <- 
            apply( cbind( object@peakChr, object@peakStart, object@peakEnd ), 1, 
                function(x) paste( x, collapse="_" ) )
        
        # summary
        
        new( "DpeakFit",
            fits=fits, optFit=optFit, optMu=optMu, optPi=optPi, optPi0=optPi0, 
            optGamma=optGamma, optDelta=optDelta, optSigma=optSigma,
            bicVec=bicVec, aicVec=aicVec, fragSet=object@fragSet, PET=object@PET, 
            fragLenTable=object@fragLenTable, Fratio=object@Fratio,
            aveFragLen=object@aveFragLen, stackedFragment=object@stackedFragment,
            peakChr=object@peakChr, peakStart=object@peakStart, peakEnd=object@peakEnd,
            estDeltaSigma=estDeltaSigma, nTop=nTop, lbDelta=lbDelta, lbSigma=lbSigma,
            psize=psize, maxComp=maxComp, pConst=pConst,
            iterInit=iterInit, iterMain=iterMain, epsilon=epsilon )
    }
)

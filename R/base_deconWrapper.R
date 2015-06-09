# deconvolution

.deconWrapper <- function( fData, estDelta=TRUE, lbDelta=c(25,25), lbSigma=c(25,25),
    psize=21, max_comp=5, pConst=0.2,
    seed=12345, niter_init=25, niter_gen=25, 
    PET, L_table=NA, Fratio=0.5, aveFragLen=NA, stop_eps=1e-6, verbose=FALSE ) {

    frag <- fData$frag
    peak <- fData$peak
    
    ##if ( is.na(frag[1,1]) ) {
    if(is.null(frag)){
    #if ( nrow(frag) == 0 ) {
        # if peak has no fragment, skip the model fitting
        
        #message( "Warning: No fragment in this peak region ",peak[1],"-",peak[2],"!" )
        success <- FALSE
    } else {           
        # starting & end positions of each peak
        ## S <- frag[,1]                   
        S <- start(frag)
        ## E <- frag[,2]
        E <- end(frag)
        
        midp <- ( S + E ) / 2
        
        if ( PET == FALSE ) {
            # SET
            
            ## strand <- frag[,3]
            strand <- as.character(strand(frag))
            ## N <- nrow(frag)
            N <- length(frag)
            
            ## locF <- which( strand=="F" )
            ## locR <- which( strand=="R" )
            locF <- which(strand == "+")
            locR <- which(strand == "-")
            FRvec <- rep( 1, N )
            FRvec[ locF ] <- Fratio
            FRvec[ locR ] <- 1 - Fratio
            
            indF <- indR <- rep( 0, N )
            if ( length(locF) > 0 ) {
                indF[ locF ] <- 1
            }
            if ( length(locR) > 0 ) {
                indR[ locR ] <- 1
            }
            
            sindex <- list( FRvec=FRvec, indF = indF, indR = indR )
            
            beta <- 2 * aveFragLen            
        } else {
            # PET
            
            strand <- NA
            L <- E - S + 1            
            fragRange <- IRanges( start=S, end=E )    
            SEL <- cbind( S, E, L )
            alpha <- .ff_dlength( L, as.numeric(names(L_table)), as.numeric(L_table) )  
        }
        
        fit_list <- mu_list <- pi_list <- pi0_list <- vector( "list", max_comp )
        
        if ( PET == FALSE ) {
            delta_list <- sigma_list <- vector( "list", max_comp )
        } else {
            gamma_list <- vector( "list", max_comp )
        }
        
        BIC_vec <- AIC_vec <- rep( NA, max_comp )
        
        for ( n_comp in 1:max_comp ) {
            # initial value: equal grid
          
            #mu_init <- seq( min((S+E)/2), max((S+E)/2), length=(n_comp+2) )
            mu_init <- seq( min(midp), max(midp), length=(n_comp+2) )
            mu_init <- mu_init[ -c(1,length(mu_init)) ]
            
            # EM algorithm
            
            set.seed( seed + n_comp )
            
            if ( PET == FALSE ) {
                # ---------------------------- SET ----------------------------
                
                # SECM
                
                fit_init <- .deconInitSET( S=S, E=E, strand=strand, peak=peak, 
                    estDelta=estDelta, lbDelta=lbDelta[1], lbSigma=lbSigma[1],
                    psize=psize, Fratio=Fratio, sindex=sindex,
                    beta=beta, niter=niter_init, mu_init=mu_init, 
                    L_table=L_table, stop_eps=stop_eps, verbose=verbose )
                
                # ECM
                
                fit <- .deconCoreSET( S=S, E=E, strand=strand, peak=peak, 
                    estDelta=estDelta, lbDelta=lbDelta[2], lbSigma=lbSigma[2],
                    psize=psize, Fratio=Fratio, sindex=sindex, 
                    beta=beta, niter=niter_gen, mu_init=fit_init$mu, 
                    pi_init=fit_init$pi, pi0_init=fit_init$pi0,
                    delta_init=fit_init$delta, sigma_init=fit_init$sigma,
                    L_table=L_table, stop_eps=stop_eps, verbose=verbose )
            } else {
                # ---------------------------- PET ----------------------------
                
                # SECM
                
                fit_init <- .deconInitPET( S=S, E=E, midp=midp, L=L, 
                    fragRange=fragRange, SEL=SEL, peak=peak, 
                    psize=psize, niter=niter_init, mu_init=mu_init, 
                    L_table=L_table, stop_eps=stop_eps, verbose=verbose )

                # ECM
                
                fit <- .deconCorePET( S=S, E=E, L=L, fragRange=fragRange, peak=peak, 
                    psize=psize, niter=niter_gen, 
                    mu_init=fit_init$mu, pi_init=fit_init$pi, pi0_init=fit_init$pi0,
                    gamma_init=fit_init$gamma,
                    L_table=L_table, alpha=alpha, stop_eps=stop_eps, verbose=verbose )
            }
            
            # keep results
            
            fit_list[[n_comp]] <- fit
            
            mu_list[[n_comp]] <- round( fit$mu )
            pi_list[[n_comp]] <- fit$pi     
            pi0_list[[n_comp]] <- fit$pi0       
            
            if ( PET == FALSE ) {  
                #pi0_list[[n_comp]] <- fit$pi0
                delta_list[[n_comp]] <- round( fit$delta )
                sigma_list[[n_comp]] <- round( fit$sigma )
            } else {
                gamma_list[[n_comp]] <- fit$gamma
            }
                 
            BIC_vec[ n_comp ] <- fit$BIC
            AIC_vec[ n_comp ] <- fit$AIC
        }
        
        success <- TRUE
    }
    
    # model selection
    
    if ( success ) {
        # if succeeded to fit the model, combine models with same dim
        
        if ( PET == FALSE ) {
            cmodel <- .combineModelSET( fit_list, mu_list, pi_list, pi0_list,
                delta_list, sigma_list, BIC_vec, AIC_vec, max_comp )
        } else {
            #cmodel <- .combineModelPET( fit_list, mu_list, pi_list, 
            #    gamma_list, BIC_vec, AIC_vec, max_comp )
            cmodel <- .combineModelPET( fit_list, mu_list, pi_list, pi0_list,
                gamma_list, BIC_vec, AIC_vec, max_comp )
        }
            
        fit_list <- cmodel$fit_list
        
        mu_list <- cmodel$mu_list
        pi_list <- cmodel$pi_list
        pi0_list <- cmodel$pi0_list
        
        if ( PET == FALSE ) {
            #pi0_list <- cmodel$pi0_list
            delta_list <- cmodel$delta_list
            sigma_list <- cmodel$sigma_list
        } else {
            gamma_list <- cmodel$gamma_list        
        }
        
        BIC_vec <- cmodel$BIC_vec
        AIC_vec <- cmodel$AIC_vec
        
        # do model selection
        
        optModel <- .selectModel( BIC_vec, max_comp, pConst )
        
        optFit <- fit_list[[ optModel ]]
        
        optMu <- mu_list[[ optModel ]]
        optPi <- pi_list[[ optModel ]]
        optPi0 <- pi0_list[[ optModel ]]
        
        if ( PET == FALSE ) {
            #optPi0 <- pi0_list[[ optModel ]]
            optDelta <- delta_list[[ optModel ]]
            optSigma <- sigma_list[[ optModel ]]
            optGamma <- NA
        } else {            
            optGamma <- gamma_list[[ optModel ]]
            #optPi0 <- optDelta <- optSigma <- NA
            optDelta <- optSigma <- NA
        }        
        
        bicVec <- BIC_vec
        aicVec <- AIC_vec
    } else {
        # if failed to fit the model, do nothing
        
        fit_list <- optFit <- optMu <- optPi <- optPi0 <-
            optGamma <- optDelta <- optSigma <- bicVec <- aicVec <- NA
    }
    
    return( list( fits=fit_list, optFit=optFit, 
        optMu=optMu, optPi=optPi, optPi0=optPi0, optGamma=optGamma, 
        optDelta=optDelta, optSigma=optSigma, bicVec=bicVec, aicVec=aicVec ) )
}

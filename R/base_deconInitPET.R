# Initial value for deconvolution, based on the heuristic EM algorithm

.deconInitPET <- function( S, E, midp, L, fragRange, SEL, peak, 
    psize=21, niter=50, mu_init, 
    L_table, stop_eps=1e-6, verbose=FALSE ) {
  
    # construct grid  
    
    grid_min <- peak[1]
    grid_max <- peak[2]
        # search binding site only within peak region
    grid_vec <- c(grid_min:grid_max)
    
    # initialization
    
    N <- length(S)
    n_group <- length(mu_init)  
        # pi for background: 0.1
    R <- grid_max - grid_min + 1
    
    mu <- mu_init
    #pi <- rep( 1/n_group, n_group )
    pi <- rep( 0.90/n_group, n_group )
    pi0 <- 0.10
    gamma <- 0.1
    
    # EM step
    
    mu_mat <- matrix( NA, niter, n_group )
    mu_mat[1,] <- mu_init
    
    pi_mat <- matrix( NA, niter, n_group )
    pi_mat[1,] <- pi
  
    pi0_vec <- rep( NA, niter )
    pi0_vec[1] <- pi0
  
    gamma_vec <- rep( NA, niter )
    gamma_vec[1] <- gamma
  
    for ( i in 2:niter ) {
        if ( verbose ) {
            print( paste("------------ iteration:",i,"------------") )
        }
    
        ########################################################################
        #                                                                      #
        #                       E step: update Z                               #
        #                                                                      #
        ########################################################################
        
        #print("E step")
        
        Z <- matrix( NA, N, n_group )        
        if ( i < niter / 2 ) {
            # normal approximation for the first half iterations
            
            for ( g in 1:n_group ) {
                Z[,g] <- pi[g] * dnorm( midp, mu[g], L/sqrt(12) )
                
                # check at least one element in Z[,g] is non-zero
                
                if ( verbose ) {
                    if ( sum(Z[,g]) == 0 ) {
                        message( "Warning: all elements in Z vector is zero!" )
                        message( "peak region: ", grid_min, "-", grid_max )
                        message( "event number: ", g )
                    }
                }
            }
        } else {
            # uniform mixture for the last half iterations
            
            muRange <- IRanges( start=mu, end=mu )
            mm <- as.matrix( findOverlaps( muRange, fragRange ) )
            mms <- split( mm[,2], mm[,1] )
            indg_list <- lapply( mms, function(x) {
                out <- rep( 0, N )
                out[x] <- 1
                return(out)
            } )
            
            for ( g in 1:n_group ) {
                if ( any( names(mms) == g ) ) { 
                    indg <- indg_list[[ as.character(g) ]]
                } else {
                    indg <- rep( 0, N )
                }
                Z[,g] <- pi[g] * ( ( ( 1 - gamma ) / L )*indg + ( gamma / (R-1) )*( 1 - indg ) )
                
                # check at least one element in Z[,g] is non-zero
                
                if ( verbose ) {
                    if ( sum(Z[,g]) == 0 ) {
                        message( "Warning: all elements in Z vector is zero!" )
                        message( "peak region: ", grid_min, "-", grid_max )
                        message( "event number: ", g )
                    }
                }
            }
        }
        
        Z0 <- pi0 / ( R + L - 1 )
        
        nMax <- .ff_ismaxZ0( Z0, Z )
        if ( nMax == N ) {   
            group <- .ff_samp( Z ) + 1
        } else {         
            group <- .ff_samp( cbind(Z0,Z) )
        }   
    
        ########################################################################
        #                                                                      #
        #                      CM step: update mu                              #
        #                                                                      #
        ########################################################################
            
        #print("M step")
        
        # group index
        
        gindex <- split( 1:N, group )
        gbinary <- lapply( gindex, function(x) {
            out <- rep( 0, N )
            out[ x ] <- 1
            return(out)
        } )
        
        glen <- lapply( gindex, length )
        gmat <- matrix( 0, n_group+1, 2 )
        gmat[,1] <- 0:n_group
        gmat[ as.numeric(names(glen))+1, 2 ] <- unlist(glen)
        
        # M step: update mu
                
        for ( g in 1:n_group ) {           
            if ( gmat[ (g+1), 2] > 0 ) {
                indg <- gindex[[ as.character(g) ]]
                SELg <- SEL[ indg, , drop=FALSE ]
                yvar <- .ff_score( grid_vec, SELg[ , 1 ], SELg[ , 2 ], SELg[ , 3 ],
                    rep( 1, gmat[(g+1),2] ), R, gamma  )
                mu_max <- grid_vec[ yvar == max(yvar) ]
                mu[g] <- median(mu_max) 
            } else {
                mu[g] <- NA
            }
        }    
            
        # M step: update pi 
        
        pi <- rep( NA, n_group )
        
        for ( g in 1:n_group ) {
            pi[g] <- gmat[ g+1, 2 ] / N 
        }
        
        pi0 <- gmat[ 1, 2 ] / N 
        
        # safe guard for pi0: when signal is weak, do not use pi0
        
        if ( pi0 > max(pi) ) {
            pi0 <- 0
            pi <- pi / sum(pi)
        }
        
        # reduce dim
        
        pi <- pi[ !is.na(mu) ]
        mu <- mu[ !is.na(mu) ]    
        n_group <- length(which( !is.na(mu) ))
        
        # safe guard for mu estimates (case: nothing left)
        
        if ( all(is.na(mu)) ) {
	        mu_old <- mu_mat[ (i-1), ]
	        mu_old <- mu_old[ !is.na(mu_old) ]
	        
	        n_group <- length(mu_old)	        
	        mu <- sample( (S+E)/2, n_group, replace=FALSE )
	        gamma <- 0.1
		    pi <- rep( 0.90/n_group, n_group )
		    pi0 <- 0.10
		    
		    next;
        }
        
        # M step: update gamma
            
        muRange <- IRanges( start=mu, end=mu )
        mm <- as.matrix( findOverlaps( muRange, fragRange ) )
        mms <- split( mm[,2], mm[,1] )
        indg_list <- lapply( mms, function(x) {
            out <- rep( 0, N )
            out[x] <- 1
            return(out)
        } )
    
        gamma <- 0
        for ( g in 1:n_group ) {
            if ( any( names(mms) == g ) ) { 
                indg <- indg_list[[ as.character(g) ]]
            } else {
                indg <- rep( 0, N )
            }
            if ( gmat[ (g+1), 2] > 0 ) {
                gamma <- gamma + sum( gbinary[[ as.character(g) ]] * ( 1 - indg ) )
            }
        }
        gamma <- gamma / N
    
        # safe guard: prevent NaN in calculating score
        
        if ( gamma < 0.01 ) {
            gamma <- 0.01
        } else if ( gamma > 0.50 ) {
            gamma <- 0.50
        }
    
        ########################################################################
        #                                                                      #
        #      Identifiability, over-fitting, track estimates & loglik         #
        #                                                                      #
        ########################################################################
            
        # identifiability problem: order constraint on \mu values
            
        pi <- pi[ order(mu) ]
        mu <- mu[ order(mu) ]
        
        # check over-fitting -> reduce dimension
        # (avoid identifiability problem due to over-fitting)
        # condition: distance <= psize
    
        if ( n_group >= 2 ) {
            mu_new <- pi_new <- c()
            mu_current <- mu[1]
            pi_current <- pi[1]
                
            for ( g in 2:n_group ) {
                if ( abs( mu[g] - mu_current ) <= psize ) {
                    mu_current <- ( mu_current + mu[g] ) / 2
                    pi_current <- pi_current + pi[g]
                } else {
                    mu_new <- c( mu_new, mu_current )
                    pi_new <- c( pi_new, pi_current )
                    
                    mu_current <- mu[g]
                    pi_current <- pi[g]
                }
            }

            mu_new <- c( mu_new, mu_current )
            pi_new <- c( pi_new, pi_current )
            
            mu <- mu_new
            pi <- pi_new
            n_group <- length(mu)
            
            # check over-fitting \pi < 0.01 -> reduce dimension
            # (avoid identifiability problem due to over-fitting)
            
            if ( any( pi < 0.01 ) ) {      
                # safeguard: reduce dim only if there is at least one remaining component
                #             if nothing remains, just stop
          
                if ( length(which( pi > 0.01 )) > 0 ) {
                    mu <- mu[ pi > 0.01 ]
                    pi <- pi[ pi > 0.01 ]      
                    n_group <- length(mu)
                } else {
      
                    # use estimates in the last iteration
                    
                    mu <- mu_mat[ (i-1), !is.na(mu_mat[(i-1),]) ]
                    mu_mat[ i, ] <- NA
                    mu_mat[ i, 1:length(mu) ] <- mu
                    
                    pi <- pi_mat[ (i-1), !is.na(pi_mat[(i-1),]) ]
                    pi_mat[ i, ] <- NA
                    pi_mat[ i, 1:length(pi) ] <- pi
            
                    pi0 <- pi0_vec[(i-1)]
                    pi0_vec[i] <- pi0
                    
                    gamma <- gamma_vec[(i-1)]
                    gamma_vec[i] <- gamma
                    
                    # stop iteration
                        
                    mu_mat <- mu_mat[ 1:i, , drop=FALSE ]
                    pi_mat <- pi_mat[ 1:i, , drop=FALSE ]
                    pi0_vec <- pi0_vec[ 1:i ]
                    gamma_vec <- gamma_vec[ 1:i ]
                    #ll <- ll[ 1:i ]
                    
                    break;
                }
            }
        
            pi <- pi / sum(pi)
        }
    
        # safeguard: if only one component & pi decrases, just stop
    
        if ( n_group == 1 & pi[1] < 0.01 ) {
            
            # use estimates in the last iteration
            
            mu <- mu_mat[ (i-1), !is.na(mu_mat[(i-1),]) ]
            mu_mat[ i, ] <- NA
            mu_mat[ i, 1:length(mu) ] <- mu
            
            pi <- pi_mat[ (i-1), !is.na(pi_mat[(i-1),]) ]
            pi_mat[ i, ] <- NA
            pi_mat[ i, 1:length(pi) ] <- pi
            
            pi0 <- pi0_vec[(i-1)]
            pi0_vec[i] <- pi0
            
            gamma <- gamma_vec[(i-1)]
            gamma_vec[i] <- gamma
            
            # stop iteration
                
            mu_mat <- mu_mat[ 1:i, , drop=FALSE ]
            pi_mat <- pi_mat[ 1:i, , drop=FALSE ]
            pi0_vec <- pi0_vec[ 1:i ]
            gamma_vec <- gamma_vec[ 1:i ]
            #ll <- ll[ 1:i ]
            break;
        
        }
        
        # track estimates
        
        mu_mat[ i, 1:length(mu) ] <- mu    
        pi_mat[ i, 1:length(pi) ] <- pi
        pi0_vec[i] <- pi0
        gamma_vec[i] <- gamma
        
        if ( verbose ) {
            print( "mu: " )
            print( mu )
            print( "pi: " )
            print( pi )
            print( "pi0: " )
            print( pi0 )
            print( "gamma: " )
            print( gamma )
        }
    } 
    
    return( list( mu=mu, pi=pi, pi0=pi0, gamma=gamma ) )
}

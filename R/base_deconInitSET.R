# Initial value for deconvolution, based on the heuristic EM algorithm

.deconInitSET <- function( S, E, strand, peak, 
    estDelta=TRUE, lbDelta=25, lbSigma=25,
    psize=21, Fratio=0.5, sindex, beta,
    niter=50, mu_init, L_table, stop_eps=1e-6, verbose=FALSE ) {
    
    # initialization
  
    grid_min <- peak[1]
    grid_max <- peak[2]
  
    N <- length(S)
    n_group <- length(mu_init)
    
    FRvec <- sindex$FRvec
    indF <- sindex$indF
    indR <- sindex$indR
        
    mu <- mu_init  
    pi <- rep( 0.90/n_group, n_group )
    pi0 <- 0.10
    
    fraglen <- E[1] - S[1] + 1
    delta <- fraglen / 2
    sigma <- delta / 2
    
    R <- grid_max - grid_min + 1    
  
    # iteration
    
    mu_mat <- matrix( NA, niter, n_group )
    mu_mat[1,] <- mu_init
    pi_mat <- matrix( NA, niter, n_group )
    pi_mat[1,] <- pi
    pi0_vec <- rep( NA, niter )
    pi0_vec[1] <- pi0
    delta_vec <- rep( NA, niter )
    delta_vec[1] <- delta
    sigma_vec <- rep( NA, niter )
    sigma_vec[1] <- sigma
  
    for ( i in 2:niter ) {    
        if ( verbose ) {
            print( paste("------------ iteration:",i,"------------") )
        }
    
        # SE step: stochastic classification
    
        Z <- matrix( 0, N, n_group )
        for ( g in 1:n_group ) {
            Z[,g] <- pi[g] * ( Fratio * dnorm( S, mu[g] - delta, sigma ) * indF +
                ( 1 - Fratio ) * dnorm( E, mu[g] + delta, sigma ) * indR )
      
            # check at least one element in Z[,g] is non-zero
            
            if ( verbose ) {
                if ( sum(Z[,g]) == 0 ) {
                    message( "Warning: all elements in Z vector is zero!" )
                    message( "peak region: ", grid_min, "-", grid_max )
                    message( "event number: ", g )
                }
            }
        }
        
        Z0 <- FRvec * pi0 / ( R + beta - 1 )
                
        nMax <- .ff_ismaxZ0( Z0, Z )
        if ( nMax == N ) {   
            group <- .ff_samp( Z ) + 1
        } else {         
            group <- .ff_samp( cbind(Z0,Z) )
        }   
        
        # CM step: update mu maximizing overlap with fragments in each group
        
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
        
        mu <- rep( NA, n_group )
        
        for ( g in 1:n_group ) {
            if ( gmat[ (g+1), 2] > 0 ) {
                indg <- gbinary[[ as.character(g) ]]
            } else {
                indg <- rep( 0, N )
            }
            
            if ( sum(indg) > 0 ) {
                sumF <- sum( ( S + delta ) * indF * indg )
                sumR <- sum( ( E - delta ) * indR * indg )
                
                mu[g] <- ( sumF + sumR ) / sum(indg)
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
		    delta <- fraglen / 2
		    sigma <- delta / 2
		    pi <- rep( 0.90/n_group, n_group )
		    pi0 <- 0.10
		    
		    next;
        }
        
        # M step: update delta
    
        if ( estDelta == TRUE ) {
            delta <- 0
            
            for ( g in 1:n_group ) {
                if ( gmat[ (g+1), 2] > 0 ) {
                    indg <- gbinary[[ as.character(g) ]]
                } else {
                    indg <- rep( 0, N )
                }
                
                if ( sum(indg) > 0 ) {
                    delta <- delta + 
                        sum( ( ( mu[g] - S ) * indF * indg ) + ( ( E - mu[g] ) * indR * indg ) )
                }
            }
            
            delta <- delta / N
        }
        
        # M step: update sigma
        
        sigma2 <- 0
        
        for ( g in 1:n_group ) {
            if ( gmat[ (g+1), 2] > 0 ) {
                indg <- gbinary[[ as.character(g) ]]
            } else {
                indg <- rep( 0, N )
            }
            
            if ( sum(indg) > 0 ) {
                sigma2 <- sigma2 + sum( ( ( S - mu[g] + delta )^2 * indF * indg ) + 
                    ( ( E - mu[g] - delta )^2 * indR *indg ) )
            }
        }
        
        sigma <- sqrt( sigma2 / N )
        
        # safe guard for delta & sigma
        
        #if ( delta < 25 ) {
        #    delta <- 25
        #}
        #if ( sigma < 25 ) {
        #    sigma <- 25
        #}
        if ( delta < lbDelta ) {
            delta <- lbDelta
        }
        if ( sigma < lbSigma ) {
            sigma <- lbSigma
        }
        
        # merge nearby mu
        
        if ( n_group >= 2 ) {
            mu_new <- c()
            pi_new <- c()
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
        }
        #print(mu)
        
        # track estimates
      
        mu_mat[ i, 1:length(mu) ] <- mu
        pi_mat[ i, 1:length(pi) ] <- pi
        pi0_vec[i] <- pi0
        delta_vec[i] <- delta
        sigma_vec[i] <- sigma
    
        if ( verbose ) {
            print( "mu: " )
            print( mu )
            print( "pi: " )
            print( pi )
            print( "pi0: " )
            print( pi0 )
            print( "delta: " )
            print( delta )
            print( "sigma: " )
            print( sigma )
        }
    }
    
    return( list( mu=mu, pi=pi, pi0=pi0, delta=delta, sigma=sigma ) )
}

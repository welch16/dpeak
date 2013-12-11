# Deconvolution, based on the generative model

.deconCoreSET <- function( S, E, strand, peak, 
    estDelta=TRUE, lbDelta=25, lbSigma=25,
    psize=21, Fratio=0.5, sindex, beta,
    niter=50, mu_init, pi_init, pi0_init, delta_init, sigma_init,
    PET=TRUE, L_table, stop_eps=1e-6, verbose=FALSE ) {  
  
    # construct grid  
    
    grid_min <- peak[1]
    grid_max <- peak[2]
        # search binding site only within peak region
    grid_vec <- c(grid_min:grid_max)
    
    # initialization
    
    N <- length(S)
    n_group <- length(mu_init)  
    L <- E - S + 1
    R <- grid_max - grid_min + 1
    
    FRvec <- sindex$FRvec
    indF <- sindex$indF
    indR <- sindex$indR
    
    mu <- mu_init
    pi <- pi_init
    pi0 <- pi0_init
    delta <- delta_init
    sigma <- sigma_init
    alpha <- rep( 1, N )  
    
    # EM step
    
    mu_mat <- matrix( NA, niter, n_group )
    mu_mat[1,] <- mu_init
    
    pi_mat <- matrix( NA, niter, n_group )
    pi_mat[1,] <- pi_init
    
    pi0_vec <- rep( NA, niter )
    pi0_vec[1] <- pi0_init
    
    delta_vec <- rep( NA, niter )
    delta_vec[1] <- delta_init
    
    sigma_vec <- rep( NA, niter )
    sigma_vec[1] <- sigma_init
    
    ll <- rep( NA, niter )
    ll[1] <- -Inf
    
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
        for ( g in 1:n_group ) {
          Z[,g] <- pi[g] * ( ( Fratio * dnorm( S, mu[g] - delta, sigma ) * indF +
                ( 1 - Fratio ) * dnorm( E, mu[g] + delta, sigma ) * indR ) )
        }
        
        Z0 <- FRvec * pi0 / ( R + beta - 1 )
        
        Znorm <- .ff_normalize( cbind(Z0,Z) )
        Z0 <- Znorm[ , 1 ]
        Z <- Znorm[ , -1, drop=FALSE ]
        
        ########################################################################
        #                                                                      #
        #                      CM step: update mu                              #
        #                                                                      #
        ########################################################################
            
        #print("M step")
        
        # M step: update mu
        
        for ( g in 1:n_group ) {
            sumF <- sum( Z[,g] * ( S + delta ) * indF )
            sumR <- sum( Z[,g] * ( E - delta ) * indR )
                    
            mu[g] <- ( sumF + sumR ) / sum( Z[,g] )
        }
            
        # M step: update pi 
        
        for ( g in 1:n_group ) {
            pi[g] <- sum( Z[,g] ) / N
        }
        pi0 <- sum( Z0 ) / N
        
        # safe guard for pi0: when signal is weak, do not use pi0
        
        if ( pi0 > max(pi) ) {
            pi0 <- 0
            pi <- pi / sum(pi)
        }
        
        # M step: update delta
    
        if ( estDelta == TRUE ) {
            delta <- 0
            
            for ( g in 1:n_group ) {
                delta <- delta + 
                    sum( Z[,g] * ( ( mu[g] - S ) * indF + ( E - mu[g] ) * indR ) )
            }
            
            delta <- delta / N
        }
        
        # M step: update sigma
        
        sigma2 <- 0
        
        for ( g in 1:n_group ) {
            sigma2 <- sigma2 + sum( Z[,g] * ( ( S - mu[g] + delta )^2 * indF + 
                ( E - mu[g] - delta )^2 * indR ) )
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
        
        if ( n_group >= 2 ) {
        
            # condition: distance <= psize
            
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
                mu <- mu[ pi > 0.01 ]
                pi <- pi[ pi > 0.01 ]      
                n_group <- length(mu)
            }
        
            norm_const <- ( 1 - pi0 ) / sum(pi)
            pi <- pi * norm_const
        }
        
        # safeguard: if only one component & pi decrases, just stop
        
        if ( n_group == 1 & pi[1] < 0.10 ) {
          
          # use estimates in the last iteration
          
          mu <- mu_mat[ (i-1), !is.na(mu_mat[(i-1),]) ]
          mu_mat[ i, ] <- NA
          mu_mat[ i, 1:length(mu) ] <- mu
          
          pi <- pi_mat[ (i-1), !is.na(pi_mat[(i-1),]) ]
          pi_mat[ i, ] <- NA
          pi_mat[ i, 1:length(pi) ] <- pi
          
          pi0 <- pi0_vec[(i-1)]
          pi0_vec[i] <- pi0
          
          delta <- delta_vec[(i-1)]
          delta_vec[i] <- delta
          
          sigma <- sigma_vec[(i-1)]
          sigma_vec[i] <- sigma
          
          # stop iteration
                
          mu_mat <- mu_mat[ 1:i, , drop=FALSE ]
          pi_mat <- pi_mat[ 1:i, , drop=FALSE ]
          pi0_vec <- pi0_vec[ 1:i ]
          delta_vec <- delta_vec[ 1:i ]
          sigma_vec <- sigma_vec[ 1:i ]
          ll <- ll[ 1:i ]
          break;
            
        }
        
        # track estimates
        
        mu_mat[ i, 1:length(mu) ] <- mu
        pi_mat[ i, 1:length(pi) ] <- pi
        pi0_vec[i] <- pi0    
        delta_vec[i] <- delta
        sigma_vec[i] <- sigma
        
        if ( verbose ) {
            print( "mu:" )
            print( mu )
            print( "pi:" )
            print( pi )
            print( "pi0:" )
            print( pi0 )
            print( "delta:" )
            print( delta )
            print( "sigma:" )
            print( sigma )
        }
        
        # track complete log likelihood
        
        ll[i] <- .loglikSET( S, E, strand, mu, pi, pi0, delta, sigma, 
            Fratio, sindex, beta, R, alpha )
        #print(ll[i])
        if ( verbose ) {
          print( "increment in loglik:" )
          print( ll[i]-ll[(i-1)] )
        }
        
        # if loglik decreases, stop iteration
        
        #print(ll[i])
        if ( ll[i] < ll[(i-1)] ) {
          #print(i)
          #print( "log lik decreases!!!" )
          
          # use estimates in the last iteration
          
          mu <- mu_mat[ (i-1), !is.na(mu_mat[(i-1),]) ]
          mu_mat[ i, ] <- NA
          mu_mat[ i, 1:length(mu) ] <- mu
          
          pi <- pi_mat[ (i-1), !is.na(pi_mat[(i-1),]) ]
          pi_mat[ i, ] <- NA
          pi_mat[ i, 1:length(pi) ] <- pi
          
          pi0 <- pi0_vec[(i-1)]
          pi0_vec[i] <- pi0
          
          delta <- delta_vec[(i-1)]
          delta_vec[i] <- delta
          
          sigma <- sigma_vec[(i-1)]
          sigma_vec[i] <- sigma
          
          # stop iteration
                
          mu_mat <- mu_mat[ 1:i, , drop=FALSE ]
          pi_mat <- pi_mat[ 1:i, , drop=FALSE ]
          pi0_vec <- pi0_vec[ 1:i ]
          delta_vec <- delta_vec[ 1:i ]
          sigma_vec <- sigma_vec[ 1:i ]
          ll <- ll[ 1:i ]
          break;
        }
        
        # check whether to stop iterations
        
        if ( ll[i] - ll[(i-1)] < stop_eps ) {
          # stop if no improvement in loglik
          
          if ( verbose ) {
            print( "stop because there is no improvements in likelihood." )
          }
                
          mu_mat <- mu_mat[ 1:i, , drop=FALSE ]
          pi_mat <- pi_mat[ 1:i, , drop=FALSE ]
          pi0_vec <- pi0_vec[ 1:i ]
          delta_vec <- delta_vec[ 1:i ]
          sigma_vec <- sigma_vec[ 1:i ]
          ll <- ll[ 1:i ]
          break;
        } else {    
          # stop if no improvement in estimates
          
          if ( length(which(!is.na(mu_mat[i,]))) == length(which(!is.na(mu_mat[(i-1),]))) &&
           all( abs( mu_mat[i,!is.na(mu_mat[i,])] - mu_mat[(i-1),!is.na(mu_mat[(i-1),])] ) < stop_eps ) &&
           all( abs( pi_mat[i,!is.na(pi_mat[i,])] - pi_mat[(i-1),!is.na(pi_mat[(i-1),])] ) < stop_eps ) &&
           abs( pi0_vec[i] - pi0_vec[(i-1)] ) < stop_eps &&
           abs( delta_vec[i] - delta_vec[(i-1)] ) < stop_eps ) {
              if ( verbose ) {
                print( "stop because there is no improvements in estimates." )
              }
              
              mu_mat <- mu_mat[ 1:i, , drop=FALSE ]
              pi_mat <- pi_mat[ 1:i, , drop=FALSE ]
              pi0_vec <- pi0_vec[ 1:i ]
              delta_vec <- delta_vec[ 1:i ]
              sigma_vec <- sigma_vec[ 1:i ]
              ll <- ll[ 1:i ]
              break;    
          }
        }
      } 
      
      aicValue <- -2 * .loglikSET( S, E, strand, mu, pi, pi0, delta, sigma, 
        Fratio, sindex, beta, R, alpha ) + (2*n_group+3) * 2
      bicValue <- -2 * .loglikSET( S, E, strand, mu, pi, pi0, delta, sigma, 
        Fratio, sindex, beta, R, alpha ) + (2*n_group+3) * log(N)
      
      return( list( mu=mu, pi=pi, pi0=pi0, delta=delta, sigma=sigma,
        mu_mat=mu_mat, pi_mat=pi_mat, pi_vec=pi0_vec, delta_vec=delta_vec, sigma_vec=sigma_vec,
        Z=Z, Z0=Z0, loglik=ll, AIC=aicValue, BIC=bicValue ) )
      #return( list( mu=mu, pi=pi, delta=delta, sigma=sigma,
      #  mu_mat=mu_mat, pi_mat=pi_mat, delta_vec=delta_vec, sigma_vec=sigma_vec,
      #  Z=Z, loglik=ll, AIC=aicValue, BIC=bicValue ) )
}

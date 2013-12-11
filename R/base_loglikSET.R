
# log likelihood

.loglikSET <- function( S, E, strand, mu, pi, pi0, 
    delta, sigma, Fratio, sindex, beta, R, alpha ) {
  
    # initialization
  
    #L <- E - S + 1
    #midp <- (S+E)/2
    #locF <- which( strand=="F" )
    #locR <- which( strand=="R" )
    #indF <- as.numeric( strand=="F" )
    #indR <- as.numeric( strand=="R" )
    
    #locF <- sindex$locF
    #locR <- sindex$locR
    FRvec <- sindex$FRvec
    indF <- sindex$indF
    indR <- sindex$indR
    
    N <- length(S)
    n_group <- length(mu)
    
    out <- alpha
  
    # likelihood: log P(L) is not incorporated (not affect loglik)
    
    out_g <- rep( 0, N )
  
    for ( g in 1:n_group ) {
        termF <- Fratio * dnorm( S, mu[g] - delta, sigma ) * indF
        termR <- ( 1 - Fratio ) * dnorm( E, mu[g] + delta, sigma ) * indR
        
        out_g <- out_g + pi[g] * ( termF + termR )
    }
    
    #termBG <- rep( pi0 / ( R + 2 * delta + 4 * sigma - 1 ), N )
    #termBG <- rep( pi0 / ( R + delta + 2 * sigma - 1 ), N )
    #termBG[ locF ] <- termBG[ locF ] * Fratio
    #termBG[ locR ] <- termBG[ locR ] * ( 1 - Fratio )
    #termBG <- FRvec * pi0 / ( R + delta + 2 * sigma - 1 )
    termBG <- FRvec * pi0 / ( R + beta - 1 )
    
    out_g <- out_g + termBG
    
    nonzero <- which( out_g > 0 )
    #out <- sum(log( out[ nonzero ] * out_g[ nonzero ] ))  
    out <- sum(log( ( out * out_g )[ nonzero ] ))  
    
    return(out)
}

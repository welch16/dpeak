
# log likelihood

#.loglikPET <- function( S, E, L, mu, pi, pi0, gamma, R, alpha ) {
.loglikPET <- function( fragRange, L, mu, pi, pi0, gamma, R, alpha ) {
  
  # initialization
  
  #N <- length(S)
  N <- length(L)
  n_group <- length(mu)
  out <- alpha
  
  # likelihood
            
    muRange <- IRanges( start=mu, end=mu )
    mm <- as.matrix( findOverlaps( muRange, fragRange ) )
    mms <- split( mm[,2], mm[,1] )
    indg_list <- lapply( mms, function(x) {
        out <- rep( 0, N )
        out[x] <- 1
        return(out)
    } )
  
  out_g <- rep( 0, N )
  for ( g in 1:n_group ) {
        #indg <- as.numeric( S <= mu[g] & mu[g] <= E )
        if ( any( names(mms) == g ) ) { 
            indg <- indg_list[[ as.character(g) ]]
        } else {
            indg <- rep( 0, N )
        }
        #out_g <- out_g + pi[g] * ( ( 1 - gamma ) / L )^indg * ( gamma / R )^( 1 - indg )
        #out_g <- out_g + pi[g] * ( ( ( 1 - gamma ) / L )*indg + ( gamma / R )*( 1 - indg ) )
        out_g <- out_g + pi[g] * ( ( ( 1 - gamma ) / L )*indg + ( gamma / (R-1) )*( 1 - indg ) )
  }
  
  termBG <- pi0 / ( R + L - 1 )
  out_g <- out_g + termBG
  
  nonzero <- which( out_g > 0 )
  #out <- sum(log( out[ nonzero ] * out_g[ nonzero ] ))  
  out <- sum(log( ( out * out_g )[ nonzero ] ))  
  
  return(out)
}

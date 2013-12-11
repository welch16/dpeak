
# calculate background by excluding combined binding sites

.bg_SET <- function( L, mu, R, n_group ) {
  denom <- R + L - 1 - L
  
  if ( n_group > 1 ) {
    for ( g in 2:n_group ) {
      denom <- denom - as.numeric( mu[g] - mu[(g-1)] <= L ) * abs( mu[g] - mu[(g-1)] ) -
        as.numeric( mu[g] - mu[(g-1)] > L ) * L
    }
  }
  
  denom <- pmax( denom, 1 ) # avoid non-positive values
  
  return( denom )
}

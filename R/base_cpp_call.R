
# Calculate score function to update mu (Rcpp)

.ff_score <- function( grid, S, E, L, Zg, R, gamma ) {
    stopifnot( length(grid) >= 1L, length(S) >= 1L, length(E) >= 1L,
        length(L) >= 1L, length(Zg) >= 1L, R > 0, gamma > 0,
        length(S) == length(E), length(E) == length(L), length(L) == length(Zg) )
    out <- .Call( "cpp_score", grid, S, E, L, Zg, R, gamma, PACKAGE="dpeak" )
    return(out)
}

# Stack fragments (Rcpp)

.ff_stack <- function( S, E, minx, maxx ) {
    stopifnot( length(S) >= 1L, length(E) >= 1L, length(S) == length(E) )
    out <- .Call( "cpp_stack", S, E, minx, maxx, PACKAGE="dpeak" )
    return(out)
}

# random sampling

.ff_samp <- function(X) {
    stopifnot( is.matrix(X), nrow(X) >= 1L, ncol(X) >= 1L )
    out <- .Call( "cpp_samp", X, PACKAGE="dpeak" )
    return(out)
}

# normalize

.ff_normalize <- function(X) {
    stopifnot( is.matrix(X), nrow(X) >= 1L, ncol(X) >= 1L )
    out <- .Call( "cpp_normalize", X, PACKAGE="dpeak" )
    return(out)
}

# match with table of fragment length

.ff_dlength <- function( L, Lname, Lfreq ) {
    stopifnot( length(L) >= 1L, length(Lname) >= 1L, length(Lfreq) >= 1L,
        length(Lname) == length(Lfreq) )
    out <- .Call( "cpp_dlength", L, Lname, Lfreq, PACKAGE="dpeak" )
    return(out)
}

#  Is Z0 larger than each Z? 

.ff_ismaxZ0 <- function( Z0, Z ) {
    stopifnot( is.matrix(Z), nrow(Z) >= 1L, ncol(Z) >= 1L, length(Z0) >= 1L, 
        nrow(Z) == length(Z0) )
    out <- .Call( "cpp_ismaxZ0", Z0, Z, PACKAGE="dpeak" )
    return(out)
}

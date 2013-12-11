\name{dpeakFit}
\alias{dpeakFit}
\alias{dpeakFit,DpeakData-method}
\title{
Fit dPeak model
}
\description{
Fit a deconvolution model.
}
\usage{
dpeakFit( object, ... )
\S4method{dpeakFit}{DpeakData}( object, estDelta=TRUE, lbDelta = 25, lbSigma = 25,
    psize=21, maxComp=5, pConst=0.2, nCore=8, verbose=FALSE, 
    seed=12345, iterInit=50, iterMain=25, epsilon=1e-6 )
}
\arguments{
    \item{object}{Object of class \code{DpeakData},
        dPeak data imported using method \code{dpeakRead}. }
    \item{estDelta}{Estimate peak-specific optimal read shift (delta) for SET data?
        Possible values are either \code{TRUE} (estimate optimal read shift for each peak)
        or \code{FALSE} (use average fragment length for all peaks).
        Default is \code{TRUE}. Not relevant when PET data is used. }
    \item{lbDelta}{Lower bound for delta parameter.}
    \item{lbSigma}{Lower bound for sigma parameter.}
    \item{psize}{Approximate size of the binding protein of interest. }
    \item{maxComp}{Maximum possible number of binding events in each peak region. }
    \item{pConst}{Value to determine the plateau in the BIC curve.
        Should be a value larger than zero and smaller than one. }
    \item{nCore}{Number of CPUs to be used when parallel computing is utilized. }
    \item{verbose}{Use verbose mode?
        Possible values are either \code{TRUE} (use) or \code{FALSE} (do not use). }
    \item{seed}{Random seed for initialization of parameters. }
    \item{iterInit}{Iteration number for initial estimation of binding sites. }
    \item{iterMain}{Iteration number for main estimation of binding sites. }
    \item{epsilon}{Criterion to stop iteration for binding site estimation. }
    \item{...}{ Other parameters to be passed through to generic \code{dpeakFit}.}  
}
\details{
Parallel computing can be utilized for faster computation 
if \code{parallel} package is installed.
Users can change the number of CPUs to be used by changing the argument \code{nCore}.
}
\value{
Construct \code{DpeakFit} class object.
}
\author{ Dongjun Chung }
\seealso{
\code{\link{dpeakRead}}, \code{\linkS4class{DpeakFit}}.
}
\examples{
\dontrun{
fit <- dpeakFit( data )
}
}
\keyword{models}
\keyword{methods}

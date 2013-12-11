\name{printEmpty}
\alias{printEmpty}
\alias{printEmpty,DpeakData-method}
\title{
Return the peak regions without any reads.
}
\description{
Return the data frame of the peak regions without any reads.
}
\usage{
printEmpty( object, ... )
\S4method{printEmpty}{DpeakData}( object, ... )
}
\arguments{
    \item{object}{Object of class \code{dpeakData},
        dPeak data obtained using the method \code{dpeakRead}. } 
    \item{...}{ Other parameters to be passed through to generic \code{printEmpty}.}
}
\author{ Dongjun Chung }
\seealso{
\code{\link{dpeakRead}}, \code{\linkS4class{DpeakData}}.
}
\examples{
\dontrun{
dataPET <- dpeakRead( peakfile="examplePeak.txt", readfile="examplePETRead.txt", 
    fileFormat="eland_result", PET=TRUE )
printEmpty( dataPET )
}
}
\keyword{models}
\keyword{methods}

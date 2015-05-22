\name{dpeak-package}
\alias{dpeak-package}
\alias{dpeak}
\docType{package}
\title{
dPeak (Deconvolution of Peaks in ChIP-seq Analysis)
}
\description{
This package provides functions for fitting dPeak, 
a statistical framework to deconvolve ChIP-seq peaks.
}
\details{
\tabular{ll}{
Package: \tab dpeak\cr
Type: \tab Package\cr
Version: \tab 2.0.1\cr
Date: \tab 2014-09-15\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains two main classes, \code{DpeakData} and \code{DpeakFit},
which represent dPeak data and deconvolution model fit, respectively.
This package contains two main methods, \code{dpeakRead} and \code{dpeakFit}.
\code{dpeakRead} method imports peak list and aligned read file and 
construct \code{DpeakData} class object.
\code{dpeakFit} method fits deconvolution model using \code{DpeakData} class object and
 constructs \code{DpeakFit} class object.
\code{DpeakFit} class object can be exported as text files
and can be used for the downstream analysis.
}
\author{
Dongjun Chung

Maintainer: Dongjun Chung <chungdon@stat.wisc.edu>
}
\keyword{ package }
\seealso{
\code{\link{dpeakRead}}, \code{\link{dpeakFit}}, \code{\link{export}}, 
\code{\linkS4class{DpeakData}}, \code{\linkS4class{DpeakFit}}.
}
\examples{
\dontrun{
# work flow for PET data

dataPET <- dpeakRead( peakfile="examplePeak.txt", readfile="examplePETRead.txt", 
    fileFormat="eland_result", PET=TRUE )
dataPET
plot( dataPET, filename="exPETplot.pdf" )

fitPET <- dpeakFit( dataPET )
fitPET
plot( fitPET, filename="exPETResult.pdf", plotType="fit" )
plot( fitPET, filename="exPETGOF.pdf", plotType="GOF" )

# work flow for SET data

dataSET <- dpeakRead( peakfile="examplePeak.txt", readfile="exampleSETRead.txt", 
    fileFormat="eland_result", PET=FALSE, fragLen=150 ) 
dataSET
plot( dataSET, filename="exSETplot_combined.pdf", strand=FALSE )
plot( dataSET, filename="exSETplot_strand_1.pdf", 
    strand=TRUE, extension=1, smoothing=TRUE )
plot( dataSET, filename="exSETplot_strand_150.pdf", 
    strand=TRUE, extension=150, smoothing=FALSE )

fitSET <- dpeakFit( dataSET )
fitSET
plot( fitSET, filename="exSETResult_combined.pdf", 
    plotType="fit", strand=FALSE )
plot( fitSET, filename="exSETResult_strand_1.pdf", 
    plotType="fit", strand=TRUE, extension=1, smoothing=TRUE )
plot( fitSET, filename="exSETResult_strand_150.pdf", 
    plotType="fit", strand=TRUE, extension=150, smoothing=FALSE )
plot( fitSET, filename="exSETGOF.pdf", plotType="GOF" )

# (common for both PET and SET data)

export( fitSET, type="txt", filename="result.txt" )
export( fitSET, type="bed", filename="result.bed" )
export( fitSET, type="gff", filename="result.gff" )
}
}

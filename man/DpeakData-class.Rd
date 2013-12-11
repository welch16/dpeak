\name{DpeakData-class}
\Rdversion{1.1}
\docType{class}
\alias{DpeakData-class}
\alias{plot,DpeakData,missing-method}
\alias{print,DpeakData-method}
\alias{show,DpeakData-method}
\title{Class "DpeakData" }
\description{
This class represents dPeak data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DpeakData", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fragSet}:}{Object of class \code{"list"}, 
        representing list of fragments for each peak. }
    \item{\code{PET}:}{Object of class \code{"logical"},
        representing whether it is paired-end tag (PET) or single-end tag (SET) data. }
    \item{\code{fragLenTable}:}{Object of class \code{"table"},
        representing distribution of fragment length when \code{PET=TRUE}. }
    \item{\code{aveFragLen}:}{Object of class \code{"numeric"},
        representing average fragment length when \code{PET=FALSE}. }
    \item{\code{Fratio}:}{Object of class \code{"numeric"},
        representing proportion of forward reads when \code{PET=FALSE}. }
    \item{\code{stackedFragment}:}{Object of class \code{"list"},
        representing number of fragments aligning to each genomic position. }
    \item{\code{peakChr}:}{Object of class \code{"character"},
        representing a vector of chromosome of each peak. }
    \item{\code{peakStart}:}{Object of class \code{"numeric"},
        representing a vector of start position of each peak. }
    \item{\code{peakEnd}:}{Object of class \code{"numeric"},
        representing a vector of end position of each peak. }
    \item{\code{emptyList}:}{Object of class \code{"character"},
        representing a vector of peak regions without reads. }
  }
}
\section{Methods}{
  \describe{
    \item{dpeakFit}{\code{signature(object = "DpeakData")}: fit the deconvolution model. }
    \item{plot}{\code{signature(x = "BinData", y = "missing", filename=NULL, 
        strand=FALSE, extension=1, smoothing=FALSE )}: 
        provide exploratory plots of fragments or reads in each peak region.
        Plots are exported to a PDF file (its file name is specified in \code{filename}).
        Options \code{strand}, \code{extension}, and \code{smoothing} are supported only for SET data.
        If \code{strand=TRUE}, reads are plotted in a strand-specific manner,
        where reads are extended to \code{extension} from its 5' end.
        If \code{smoothing=TRUE}, a smoothed plot (using the smoothing spline) is provided.
        If \code{strand=FALSE}, strand information is ignored.}
    \item{print}{\code{signature(x = "DpeakData")}: (not supported yet) }
    \item{printEmpty}{\code{signature(x = "DpeakData")}: 
        provide the data frame of peak regions without reads. }
    \item{show}{\code{signature(object = "DpeakData")}: provide brief summary of the object. }
     }
}
\author{ Dongjun Chung}
\seealso{
\code{\link{dpeakRead}}, \code{\link{dpeakFit}}.
}
\examples{
showClass("DpeakData")
\dontrun{
# PET data

dataPET <- dpeakRead( peakfile="examplePeak.txt", readfile="examplePETRead.txt", 
    fileFormat="eland_result", PET=TRUE )
dataPET
plot( dataPET, filename="exPETplot.pdf" )
fitPET <- dpeakFit( dataPET )

# SET data

dataSET <- dpeakRead( peakfile="examplePeak.txt", readfile="exampleSETRead.txt", 
    fileFormat="eland_result", PET=FALSE, fragLen=150 )
dataSET
plot( dataSET, filename="exSETplot_combined.pdf", strand=FALSE )
plot( dataSET, filename="exSETplot_strand_1.pdf", 
    strand=TRUE, extension=1, smoothing=TRUE )
plot( dataSET, filename="exSETplot_strand_150.pdf", 
    strand=TRUE, extension=150, smoothing=FALSE )
fitSET <- dpeakFit( dataSET )
}
}
\keyword{classes}

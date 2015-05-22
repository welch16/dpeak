\name{DpeakFit-class}
\Rdversion{1.1}
\docType{class}
\alias{DpeakFit-class}
\alias{plot,DpeakFit,missing-method}
\alias{print,DpeakFit-method}
\alias{show,DpeakFit-method}
\title{Class "DpeakFit" }
\description{
This class represents deconvolution model fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DpeakFit", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fits}:}{Object of class \code{"list"},
        representing list of deconvolution fits of all possible models for each peak. }
    \item{\code{optFit}:}{Object of class \code{"list"},
        representing list of fits of the optimal model for each peak. }
    \item{\code{optMu}:}{Object of class \code{"list"},
        representing list of binding sites of the optimal model for each peak. }
    \item{\code{optPi}:}{Object of class \code{"list"},
        representing list of relative strengths of the optimal model for each peak. }
    \item{\code{optPi0}:}{Object of class \code{"list"},
        representing list of background strengths of the optimal model for each peak. }
    \item{\code{optGamma}:}{Object of class \code{"list"},
        representing list of background proportion of the optimal model for each peak. }
    \item{\code{optDelta}:}{Object of class \code{"list"},
        representing list of read shift of the optimal model for each peak
        when \code{PET=FALSE}. }
    \item{\code{optSigma}:}{Object of class \code{"list"},
        representing list of read standard deviation of the optimal model for each peak
        when \code{PET=FALSE}. }
    \item{\code{bicVec}:}{Object of class \code{"list"},
        representing list of BIC values for each peak. }
    \item{\code{aicVec}:}{Object of class \code{"list"},
        representing list of AIC values for each peak. }
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
    \item{\code{estDeltaSigma}:}{Object of class \code{"character"},
        representing the approach to estimate delta and sigma parameters for SET data. }
    \item{\code{nTop}:}{Object of class \code{"numeric"},
        representing the number of candidate regions used to estimate common delta and sigma estimates. }
    \item{\code{lbDelta}:}{Object of class \code{"numeric"},
        representing a lower bound for the delta parameter. }
    \item{\code{lbSigma}:}{Object of class \code{"numeric"},
        representing a lower bound for the sigma parameter. }
    \item{\code{psize}:}{Object of class \code{"numeric"},
        representing approximate size of the binding protein of interest. }
    \item{\code{maxComp}:}{Object of class \code{"numeric"},
        representing maximum possible number of binding events in each peak region. }
    \item{\code{pConst}:}{Object of class \code{"numeric"},
        representing value to determine the plateau in the BIC curve. }
    \item{\code{iterInit}:}{Object of class \code{"numeric"},
        representing iteration number for initial estimation of binding sites. }
    \item{\code{iterMain}:}{Object of class \code{"numeric"},
        representing iteration number for main estimation of binding sites. }
    \item{\code{epsilon}:}{Object of class \code{"numeric"},
        representing criterion to stop iteration for binding site estimation. }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "DpeakFit", y = "missing", filename=NULL, 
        plotType="fit", strand=FALSE, extension=1, smoothing=FALSE,
        threshold=0.10, nsimul=10000, seed=12345, nCore=8)}: 
        draw plots of deconvolution results if \code{plotType="fit"}, 
        goodness of fit (GOF) plots if \code{plotType="GOF"},
        or plots of Bayesian information criterion (BIC) 
        and Akaike information criterion (AIC) curves if \code{plotType="BIC"}.
        Plots are exported to a PDF file (its file name is specified in \code{filename}).
        In deconvolution result plots, binding sites with strength larger than \code{threshold}
        are drawn in dark blue and other binding sites are drawn in light blue.
        When \code{plotType="fit"}, 
        options \code{strand}, \code{extension}, and \code{smoothing} are supported for SET data.
        If \code{strand=TRUE}, reads are plotted in a strand-specific manner,
        where reads are extended to \code{extension} from its 5' end.
        If \code{smoothing=TRUE}, a smoothed plot (using the smoothing spline) is provided.
        If \code{strand=FALSE}, strand information is ignored.
        For the GOF plots, \code{nsimul} fragments are simulated from the fitted model
        (\code{seed} indicates random seed; \code{nCore} CPUs are used for parallel computing).}
    \item{print}{\code{signature(x = "DpeakFit")}: (not supported yet) }
    \item{show}{\code{signature(object = "DpeakFit")}:
        provide brief summary of the object. }
     }
}
\author{ Dongjun Chung }
\seealso{
\code{\link{dpeakFit}}, \code{\link{export}}.
}
\examples{
showClass("DpeakFit")
\dontrun{
# PET data

fitPET <- dpeakFit( dataPET )
fitPET
plot( fitPET, filename="exPETResult.pdf", plotType="fit" )
plot( fitPET, filename="exPETGOF.pdf", plotType="GOF" )

# SET data

fitSET <- dpeakFit( dataSET )
fitSET
plot( fitSET, filename="exSETResult_combined.pdf", 
    plotType="fit", strand=FALSE )
plot( fitSET, filename="exSETResult_strand_1.pdf", 
    plotType="fit", strand=TRUE, extension=1, smoothing=TRUE )
plot( fitSET, filename="exSETResult_strand_150.pdf", 
    plotType="fit", strand=TRUE, extension=150, smoothing=FALSE )
plot( fitSET, filename="exSETGOF.pdf", plotType="GOF" )
plot( fitSET, filename="exSETBIC.pdf", plotType="BIC" )}
}
\keyword{classes}

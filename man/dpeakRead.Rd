\name{dpeakRead}
\alias{dpeakRead}
\title{
Import peak list and aligned read files
}
\description{
Import and process peak list and aligned read files.
}
\usage{
dpeakRead( peakfile=NULL, readfile=NULL, fileFormat="eland_result", 
    PET=FALSE, fragLen=200, parallel=FALSE, nCore=8, tempDir=NULL, perl = "perl" )
}
\arguments{
  \item{peakfile}{File name of the peak list.}
  \item{readfile}{Name of the aligned read file.}
  \item{fileFormat}{
    Format of the aligned read file to be processed.
    For single-end tag (SET) ChIP-seq data,
    \code{dpeakRead} permits the following aligned read file formats:
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        \code{"sam"} (SAM), and \code{"bed"} (BED).
    For paired-end tag (PET) ChIP-seq data,
    \code{dpeakRead} permits only \code{"eland_result"} (Eland result format).
}
  \item{PET}{
    Is the aligned read file paired-end tag (PET)? 
    Possible values are either \code{TRUE} (PET) or \code{FALSE} (SET).
    Default is \code{FALSE} (SET).
}
  \item{fragLen}{
    Average fragment length. Default is 200. Not relevant when \code{PET=TRUE}.
}
  \item{parallel}{Utilize multiple CPUs for parallel computing 
    using \code{"parallel"} package?
    Possible values are \code{TRUE} (use \code{"parallel"}) 
    or \code{FALSE} (not use \code{"parallel"}).
    Default is \code{FALSE} (not use \code{"parallel"}).
}
  \item{nCore}{Number of CPUs when parallel computing is utilized.
}
  \item{tempDir}{
    Directory to store temporary files. If \code{tempDir=NULL}, 
    \code{dpeakRead()} will use the temporary directory used by R.
}
  \item{perl}{
    Name of the perl executable to be called. Default is \code{"perl"}.
}
}
\details{
The first three columns of the peak list file (specifed as \code{peakfile})
are assumed to be chromosome, start and end positions of each peak region.
There should be no header in the peak list file.

When the data contains multiple chromosomes, 
parallel computing can be utilized for faster preprocessing
if \code{parallel=TRUE} and \code{parallel} package is installed.
\code{nCore} determines number of CPUs used for parallel computing.
}
\value{
Construct \code{DpeakData} class object.
}
\author{ Dongjun Chung }
\seealso{
\code{\link{dpeakFit}}, \code{\linkS4class{DpeakData}}.
}
\examples{
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
fitSET <- dpeakFit( dataSET )}
}
\keyword{models}
\keyword{methods}

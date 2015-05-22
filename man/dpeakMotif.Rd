\name{dpeakMotif}
\alias{dpeakMotif}
\title{
Implement de novo motif analysis based on the peak list
}
\description{
Implement de novo motif analysis based on the peak list, using MEME and FIMO.
}
\usage{
dpeakMotif( peakfile=NULL, refGenome=NULL, flanking=100,
	memeArgument="-dna -mod zoops -nmotifs 1 -minw 10 -maxw 20 -revcomp -maxsize 1000000000", 
	fimoArgument="-max-stored-scores 100000000 -motif-pseudo 0.000001",
    tempDir=NULL )
}
\arguments{
  \item{peakfile}{File name of the peak list.}
  \item{refGenome}{\code{BSgenome} class object to extract sequences.}
  \item{flanking}{Flanking length.}
  \item{memeArgument}{Parameters for \code{MEME}.}
  \item{fimoArgument}{Parameters for \code{FIMO}.}
  \item{tempDir}{Directory of temporary files for sequence extraction, MEME, and FIMO.}
}
\details{
The first three columns of the peak list file (specifed as \code{peakfile})
are assumed to be chromosome, start and end positions of each peak region.
There should be no header in the peak list file.

\code{refGenome} is a \code{BSgenome} class object and assumed to already be available in the R environment.
}
\value{
Construct \code{DpeakMotif} class object.
}
\author{ Dongjun Chung }
\seealso{
\code{\link{dpeakFit}}, \code{\linkS4class{DpeakMotif}}.
}
\examples{
\dontrun{
library(BSgenome.Ecoli.NCBI.20080805)
resultMotif <- dpeakMotif( peakfile="examplePeak.txt", refGenome=Ecoli )
resultMotif
fitPET <- dpeakFit( dataPET, resultMotif )
}
}
\keyword{models}
\keyword{methods}


# export all possible info

#.exportTXT <- function( chrVec, muVec, piVec, psize, nameVec, filename )
.exportTXT <- function( chrVec, muVec, strengthVec, psize, nameVec, filename )
{
    outFormat <- data.frame( chrVec, round(muVec)-floor(psize/2), round(muVec)+floor(psize/2),
    #    piVec/1000, nameVec, stringsAsFactors=FALSE )
        strengthVec, nameVec, stringsAsFactors=FALSE )
    colnames(outFormat) <- c( "chrID", "siteStart", "siteEnd", "strength", "peakName" )
    
    # Note (ver 0.9.9): chromStart 1-base & chromEnd inclusive.
    outFormat[ outFormat[,2] == 0, 2 ] <- 1
    	# first base should be 1, not 0, if we use 1-base system
        
    # variable names
    
    cat( file=filename )
    cat( as.character(colnames(outFormat)), file=filename, sep="\t", append=TRUE )
    cat( "\n", file=filename, append=TRUE )
    
    # peak list
     
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
        cat( "\n", file=filename, append=TRUE )
    }
    
    message( "Info: Dpeak results were exported in TXT format:" )
    message( "Info: File name: ", filename )
}

# GFF format (score = 1000 * \pi)

#.exportGFF <- function( chrVec, muVec, piVec, psize, nameVec, filename )
.exportGFF <- function( chrVec, muVec, strengthVec, psize, nameVec, filename )
{
    # GFF: seqname, source, feature, start, end, score, strand, frame, group
    # Note (ver 0.9.9): chromStart 1-base & chromEnd inclusive.
    
    outFormat <- data.frame( chrVec, "Dpeak", nameVec, 
    #    round(muVec)-floor(psize/2), round(muVec)+floor(psize/2), piVec,
        round(muVec)-floor(psize/2), round(muVec)+floor(psize/2), strengthVec,
        ".", ".", ".", stringsAsFactors=FALSE )
    outFormat[ outFormat[,4] <= 0, 4 ] <- 1
    	# first base should be 1, not 0, if we use 1-base system
    
    line0 <- 'track name=dpeak description=\"Dpeak binding sites\"'
    cat( as.character(line0), "\n", file=filename )
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
        cat( "\n", file=filename, append=TRUE )
    }
    
    message( "Info: Dpeak results were exported in GFF format:" )
    message( "Info: File name: ", filename )
}

# BED format (score = 1000 * \pi)

#.exportBED <- function( chrVec, muVec, piVec, psize, nameVec, filename )
.exportBED <- function( chrVec, muVec, strengthVec, psize, nameVec, filename )
{
    # BED: (required) chrom, chromStart, chromEnd
    # BED: (optional) name, score, strand, thickStart, thinkEnd, itemRgb,
    #                   blockCount, blockSizes, blockStarts
    # Note (ver 0.9.9): chromStart 0-base & chromEnd NOT inclusive (i.e., 1-base).
    
    outFormat <- data.frame( chrVec, round(muVec)-floor(psize/2), round(muVec)+floor(psize/2),
    #    nameVec, piVec, stringsAsFactors=FALSE )
        nameVec, strengthVec, stringsAsFactors=FALSE )
    outFormat[ , 2 ] <- outFormat[ , 2 ] - 1
    outFormat[ outFormat[,2] < 0, 2 ] <- 0
    	# first base becomes -1 by (peakList$peakStart-1)
    
    line0 <- 'track name=dpeak description=\"Dpeak binding sites\" useScore=1'
    cat( as.character(line0), "\n", file=filename )
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
        cat( "\n", file=filename, append=TRUE )
    }
    
    message( "Info: Dpeak results were exported in BED format:" )
    message( "Info: File name: ", filename )
}

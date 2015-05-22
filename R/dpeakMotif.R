
# read bin-level data & process it

dpeakMotif <- function( peakfile=NULL, refGenome=NULL, flanking=100,
	memeArgument="-dna -mod zoops -nmotifs 1 -minw 10 -maxw 20 -revcomp -maxsize 1000000000", 
	fimoArgument="-max-stored-scores 100000000 -motif-pseudo 0.000001",
    tempDir=NULL )
{         
    # check refGenome
	
	if ( is.null(refGenome) | !is( refGenome, "BSgenome" ) ) {
		stop( "Please provide the appropriate reference genome BSgenome object!" )
	}	
	
	# check whether meme exists
	
	#CMD <- "meme 2>&1"
    #suppressWarnings( res <- system( CMD, intern = TRUE ) )
  
    #if ( length( grep( "USAGE", res ) ) == 0 ) {
        # cannot proceed if meme does not exist
        
    #    stop( "meme is not found on your system! Either check $PATH if installed or please install meme." )
    #}
	
	# check whether fimo exists
    
    #CMD <- "fimo 2>&1"
    #suppressWarnings( res <- system( CMD, intern = TRUE ) )
  
    #if ( length( grep( "USAGE", res ) ) == 0 ) {
        # cannot proceed if fimo does not exist
        
    #    stop( "fimo is not found on your system! Either check $PATH if installed or please install fimo." )
    #}
	
	# process peak set & reads (assume BED file format)
    # - peak set: assume that first 3 columns of both files are chr, start, end
    
    message( "Info: Reading peak list..." )
    peakSet <- read.table( peakfile, header=FALSE, stringsAsFactors=FALSE, sep="\t" )
    nPeak <- nrow(peakSet)
    gc()
    #print( Sys.time() )
	
	# set up temporary files	
    
    if ( is.null(tempDir) ) {
        tempDir <- tempdir()
    }
	
	tempFASTA <- paste(tempDir,"/FASTA",sep="")
	tempMEME <- paste(tempDir,"/MEME",sep="")
	tempFIMO <- paste(tempDir,"/FIMO",sep="")
	tempFimoLog <- paste(tempDir,"/fimo_log.txt",sep="")
	
	# set peak coordinates
	
	peakChr <- peakSet[,1]
	peakStart <- peakSet[,2] - flanking + 1
	peakEnd <- peakSet[,3] + flanking - 1
	
	# adjust start coordinates
	
	peakStart[ peakStart <= 0 ] <- 1
		# coordinate cannot be negative
	
	# adjust end coordinates
	
	chrlist <- unique(peakSet[,1])
	locChr <- split( 1:nrow(peakSet), peakSet[,1] )
	
	for ( i in 1:length(chrlist) ) {
		
		chr.i <- chrlist[i]
		
		# extract coordinates for each chromosome		
		
		loc.i <- locChr[[ chr.i ]]		
		#peakStart.i <- peakStart[ loc.i ]
		peakEnd.i <- peakEnd[ loc.i ]
		
		# length of chromosome
		
		len.i <- length(refGenome[[ chr.i ]])
		
		# adjust end coordinates & return
		
		peakEnd.i[ peakEnd.i > len.i ] <- len.i
		peakEnd[ loc.i ] <- peakEnd.i
	}
	
	# extract sequences corresponding to peak regions
	
	message( "Info: Extracting sequences..." )
	
	seqSet <- getSeq( x = refGenome, 
		names = peakSet[,1], start = peakStart, end = peakEnd, 
		strand = "+" )
	gc()
	
	if ( nrow(peakSet) > 1 ) {
		for ( i in 1:length(seqSet) ) {
			# ID
			
			if ( i == 1 ) {
				cat( ">seq",i,"\n", sep="", file=tempFASTA, append=FALSE )
			} else {
				cat( ">seq",i,"\n", sep="", file=tempFASTA, append=TRUE )
			}
			
			# sequence
			
			cat( as.character(seqSet[i]),"\n", sep="", file=tempFASTA, append=TRUE )
		}
	} else {
		# ID
		
		cat( ">seq",i,"\n", sep="", file=tempFASTA, append=FALSE )
		
		# sequence
		
		cat( as.character(seqSet),"\n", sep="", file=tempFASTA, append=TRUE )		
	}
	
	rm(seqSet)
	gc()
	
	# run MEME to identify de novo motifs
	
	message( "Info: Running MEME..." )
	
	CMD <- paste( "meme ",tempFASTA," -oc ",tempMEME," ",memeArgument, sep="" )
	#system( CMD, ignore.stdout=FALSE, ignore.stderr=FALSE )
	system(CMD)
	
	motifFile <- paste( tempMEME,"/meme.txt", sep="" )
	outMEME <- scan( file=motifFile, what="character", sep="\n" )
	motifIden <- outMEME[ grep( "regular expression", outMEME ) + 2 ]
	
	# run FIMO to identify locations of motifs on peak regions
	
	message( "Info: Running FIMO..." )
	
	fimoLogFile <- tempFimoLog
	
	CMD <- paste( "fimo ",fimoArgument," -oc ",tempFIMO," ",motifFile," ",tempFASTA," > ",fimoLogFile," 2>&1", sep="" )
	#system( CMD, ignore.stdout=FALSE, ignore.stderr=FALSE )
	system(CMD)
	
	# summarize FIMO results as coordinates within peaks
	
	outFIMO <- read.table( file=paste( tempFIMO,"/fimo.txt", sep="" ),
		header=TRUE, sep="\t", comment.char="", stringsAsFactors=FALSE )
		# dir/fimo.txt contains output & header starts with "#"		
	
	if ( nrow(outFIMO) > 0 ) {
		peakMotifOrg <- as.numeric( gsub( "seq", "", outFIMO$sequence.name ) )
		locMotifOrg <- round( ( outFIMO$start + outFIMO$stop ) / 2 )
		locMotifVec <- peakStart[ peakMotifOrg ] + locMotifOrg - 1	
		listMotifVec <- split( locMotifVec, peakMotifOrg )
	
		locMotif <- vector( "list", nrow(peakSet) )
		for ( i in 1:nrow(peakSet) ) {
			locMotif[[i]] <- NA
		}	
		for ( i in 1:length(listMotifVec) ) {
			locMotif[[ as.numeric(names(listMotifVec)[i]) ]] <- listMotifVec[[i]]
		}
	
		npeakWithMotif <- length(listMotifVec)
		medWidth <- median( sapply( listMotifVec, length ) )

		# info about preprocessing
		# - identified motif, # detected motifs
    
	    cat( "------------------------------------------------------------\n" )
	    cat( "Info: Preprocessing summary\n" )
	    cat( "------------------------------------------------------------\n" )
	    cat( "Identified motif:\n", sep="" )
		for ( i in 1:length(motifIden) ) {
			cat( "\t",motifIden[i],"\n", sep="" )
		}
	    cat( "Number of peaks containing detected motifs: ",npeakWithMotif,"\n", sep="" )
	    cat( "Number of peaks missing detected motifs: ",( nrow(peakSet) - npeakWithMotif ),"\n", sep="" )
	    cat( "Median number of regulatory sequences in each peak: ",medWidth,"\n", sep="" )
		cat( "\t(after excluding peaks missing detected motifs)\n" )
	    cat( "------------------------------------------------------------\n" )
	} else {
		# if there is no peak including motif
		
		locMotif <- vector( "list", nrow(peakSet) )
		for ( i in 1:nrow(peakSet) ) {
			locMotif[[i]] <- NA
		}			

		# info about preprocessing
		# - identified motif, # detected motifs
    
	    cat( "------------------------------------------------------------\n" )
	    cat( "Info: Preprocessing summary\n" )
	    cat( "------------------------------------------------------------\n" )
	    cat( "Identified motif:\n", sep="" )
		for ( i in 1:length(motifIden) ) {
			cat( "\t",motifIden[i],"\n", sep="" )
		}
		cat( "Note: No peak contains detected motifs.\n")
	    cat( "------------------------------------------------------------\n" )
	}
    
    new( "DpeakMotif", motif=motifIden, locMotif=locMotif,
        peakChr=peakChr, peakStart=peakStart, peakEnd=peakEnd )
}

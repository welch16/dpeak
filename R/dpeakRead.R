
# read bin-level data & process it

dpeakRead <- function( peakfile=NULL, readfile=NULL,
    fileFormat="eland_result", PET=FALSE, fragLen=200, 
    parallel=FALSE, nCore=8, tempDir=NULL, perl="perl" )
{            
    # process aligned read file

    ## browser()
    
    if ( PET ) {
        message( "Info: Paired-end tag (PET) is assumed (PET=TRUE)." )
    } else {
        message( "Info: Single-end tag (SET) is assumed (PET=FALSE)." )
        message( "Info: Average fragment length is set as ",fragLen," (fragLen=",fragLen,")." )
    }
    #print( Sys.time() )
    
    message( "Info: Reading and processing aligned read file..." )

    
    ## define intermediate file names (not sure if we need them though)
    if ( is.null(tempDir) ) {
        tempfileName <- tempfile( c("output","summary") )
    } else {      
        tempfileName <- c( paste(tempDir,"output.txt",sep=""),
            paste(tempDir,"summary.txt",sep="") )
    }

    readData <- .loadReadData(object = peakfile,readfile =readfile,
                  fileFormat = fileFormat,PET = PET,fragLen = fragLen,
                  keepReads = TRUE,parallel = parallel,nCore = nCore,
                  tempDir = tempDir,perl = perl)

    fragSet <- readData$fragSet
    stackedFragment <- readData$stackedFragment
        
        # intermediate file name
    ## .constructExtRead( infile=readfile, outfile=tempfileName[1],
    ##     summaryfile=tempfileName[2], fileFormat=fileFormat, 
    ##     PET=PET, fragLen=fragLen, perl=perl )            
    
    # read summary file (chrID, # lines)
    
    ## summaryInfo <- read.table( tempfileName[2], header=FALSE,
    ##     stringsAsFactors=FALSE, comment.char="", check.names=FALSE )
    ## colnames(summaryInfo) <- c("chrID","nline")   
    ## readChr <- summaryInfo[,1] 
    ## #print( Sys.time() )
        
    # process peak set & reads (assume BED file format)
    # - peak set: assume that first 3 columns of both files are chr, start, end
    
    message( "Info: Reading peak list..." )
    if( is(peakfile, "GRanges") ){
        peakSet <- as(peakfile, "data.frame")
        peakSet <- within( peakSet, seqnames <- as(seqnames, "character"))
      }else{
        peakSet <- read.table( peakfile, header=FALSE, stringsAsFactors=FALSE )
      }
    nPeak <- nrow(peakSet)
    gc()
    #print( Sys.time() )
    
    # match chromosomes between peak list and reads
    
    peakByChr <- split( peakSet[ , 2:3, drop=FALSE ], peakSet[,1] )
    peakChr <- names(peakByChr)


    peakChr <- peakSet[,1]
    peakStart <- peakSet[,2]
    peakEnd <- peakSet[,3]
    ## peakStart <- do.call(c,lapply(peakByChr,function(x)x[,1]))
    ## peakEnd <- do.call(c,lapply(peakByChr,function(x)x[,2]))

    
    ## chrCommon <- Reduce( intersect, list( peakChr, readChr ) )
    ## chrCommon <- sort(chrCommon)
    
    # match reads to each peak region
    # (using parallel computing, if parallel exists)
    
    ## message( "Info: Processing and combining peak list and reads..." )    
    
    ## if ( parallel == TRUE ) {        
    ##     out <- mclapply( chrCommon, 
    ##         function(x) .matchPeakRead( chr=x, peakCur=peakByChr[[x]],
    ##             outfileName=paste(tempfileName[1],"_",x,sep=""), 
    ##             nRow=summaryInfo[ summaryInfo[,1]==x, 2 ], PET=PET ), 
    ##         mc.cores = nCore )
    ## } else {        
    ##     out <- lapply( chrCommon, 
    ##         function(x) .matchPeakRead( chr=x, peakCur=peakByChr[[x]],
    ##             outfileName=paste(tempfileName[1],"_",x,sep=""), 
    ##             nRow=summaryInfo[ summaryInfo[,1]==x, 2 ], PET=PET )
    ##         )
    ## }   
    
    ## fragSet <- vector( "list", nPeak )
    ## peakChr <- peakStart <- peakEnd <- rep( NA, nPeak )
    ## nameVec <- rep( NA, nPeak )
    
    ## nEmpty <- 0
    ## emptyList <- c()
    ## cur <- 1
    
    ## if ( PET == FALSE ) {
    ##     nF <- nAll <- rep( NA, length(chrCommon) ) 
    ## }
    
    ## for ( chr in 1:length(chrCommon) ) {
    ##     nPeakCur <- length( out[[chr]]$nameVecCur )
        
    ##     for ( j in 1:nPeakCur ) {
    ##         if ( !is.na(out[[chr]]$fragSetCur[[j]][1,1]) ) {
    ##             fragSet[[cur]] <- out[[chr]]$fragSetCur[[j]]
    ##         } else {
    ##             fragSet[[cur]] <- matrix( NA )
    ##             nEmpty <- nEmpty + 1
    ##             emptyList <- c( emptyList, out[[chr]]$nameVecCur[j] )
    ##         }
    ##         peakChr[cur] <- out[[chr]]$peakChrCur[j]
    ##         peakStart[cur] <- out[[chr]]$peakStartCur[j]
    ##         peakEnd[cur] <- out[[chr]]$peakEndCur[j]
    ##         nameVec[cur] <- out[[chr]]$nameVecCur[j]            
    ##         cur <- cur + 1
    ##     }
        
    ##     if ( PET == FALSE ) {
    ##         nF[chr] <- out[[chr]]$nFCur
    ##         nAll[chr] <- out[[chr]]$nAllCur
    ##     }
    ## }    
    ## names(fragSet) <- nameVec
    
    ## if ( PET == TRUE ) {    
    ##     fragLen <- c()        
    ##     for ( chr in 1:length(chrCommon) ) {
    ##         fragLen <- c( fragLen, out[[chr]]$fragLenCur )
    ##     }   
    ## }

    emptyPeakIdx <- readData$numReads == 0
    nEmpty <- sum(emptyPeakIdx)

    
    if ( nEmpty == 0 ) {
        emptyList <- ""
    } else {
        emptyList <- names(which(emptyPeakIdx))
    }
    
    ## rm( out )
    ## gc()    
    ## #print( Sys.time() )
        
    # check proportion of forward reads for SET data    
    
    if ( PET == TRUE ) {
        Fratio <- 0.5
    } else {
        nAll <- sum(readData$numReads)
        nF <- mclapply(fragSet[-emptyPeakIdx],function(x)sum(as.character(strand(x))=="+"),
           mc.preschedule = TRUE,mc.cores = nCore)
        nF <- sum(unlist(nF))
        Fratio <- sum(nF) / sum(nAll)
    }
    
    # stack fragment (projection to coordinates)
    
    ## if ( parallel == TRUE ) {        
    ##     stackedFragment <- mclapply( fragSet, .stackFragment, mc.cores = nCore )
    ## } else {        
    ##     stackedFragment <- lapply( fragSet, .stackFragment )
    ## }   
    ## names(stackedFragment) <- nameVec
    #print( Sys.time() )
    
    # if PET, calculate distribution of fragment length

    
    if ( PET ) {
        message( "Info: Calculating distribution of fragment length..." )
        
        fragLenTable <- table( fragLen )
        aveFragLen <- median( fragLen )
    } else {    
        fragLenTable <- table( fragLen )
        aveFragLen <- fragLen
    }
    gc()
    
    message( "Info: Done!\n" )
    #print( Sys.time() )
    
    # remove temporary files after use
    
    unlink( tempfileName[1] )
    unlink( tempfileName[2] )
    
    # info about preprocessing
    
    nFrag <- readData$numReads

    sumRead <- sum(nFrag)
    medNumRead <- median(nFrag)
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: Preprocessing summary\n" )
    cat( "------------------------------------------------------------\n" )
    cat( "Tag type: ",ifelse(PET,"PET","SET"),"\n", sep="" )
    cat( "Number of chromosomes: ",length(peakChr),"\n", sep="" )
    cat( "Number of peaks: ",nPeak,"\n", sep="" )
    cat( "Number of peaks without reads: ",nEmpty,"\n", sep="" )
    cat( "[Note] Use 'printEmpty' method to check the list.\n" )
    cat( "Number of utilized reads: ",sumRead,"\n", sep="" )
    cat( "Median number of reads in each peak: ",medNumRead,"\n", sep="" )
    if ( PET == TRUE ) {
        cat( "Median fragment length: ",aveFragLen,"\n", sep="" )    
    } else {
        cat( "Fragment length (provided by user): ",aveFragLen,"\n", sep="" )  
        Fper <- round( 100 * Fratio )
        Rper <- 100 - Fper 
        cat( "Percentage of forward strand reads: ",Fper," %\n", sep="" ) 
        cat( "Percentage of reverse strand reads: ",Rper," %\n", sep="" )    
    }
    cat( "------------------------------------------------------------\n" )

    new( "DpeakData", fragSet=fragSet, PET=PET, fragLenTable=fragLenTable,
        aveFragLen=aveFragLen, Fratio=Fratio, stackedFragment=stackedFragment,
        peakChr=peakChr, peakStart=peakStart, peakEnd=peakEnd,
        emptyList=emptyList )
}

## # match peak & reads

## .matchPeakRead <- function( chr, peakCur, outfileName, nRow, PET ) {
    
##     # read processed read file
##     # - PET read: chr, start, end
##     # - SET read: chr, position, strand, read length
##     ## browser()
  
##     if ( PET == TRUE ) {
##         # if PET, (chrID, start, end)
        
##         readCur <- read.table( outfileName, sep='\t', nrows=nRow,
##             header=FALSE, stringsAsFactors=FALSE, comment.char="",
##             colClasses=c("character","numeric","numeric"), check.names=FALSE )
##         colnames(readCur) <- c("chrID","start","end")
##     } else {    
##         # if SET, (chrID, start, end, strand)
        
##         readCur <- read.table( outfileName, sep='\t', nrows=nRow,
##             header=FALSE, stringsAsFactors=FALSE, comment.char="",
##             colClasses=c("character","numeric","numeric","character"), 
##             check.names=FALSE )
##         colnames(readCur) <- c("chrID","start","end","str")
##     }
##     gc()        
    
##     # match reads for each peak region
    
##     peakRange <- IRanges( start=peakCur[,1], end=peakCur[,2] )
##     midpt <- ( readCur[,2] + readCur[,3] ) / 2
##     readRange <- IRanges( start=midpt, end=midpt )        
##     readMatch <- as.matrix( findOverlaps( peakRange, readRange ) )   
##     matchList <- split( readMatch[,2], readMatch[,1] )
    
##     # check whether there is any peak region without reads
    
##     existInd <- rep( 0, nrow(peakCur) )
##     existInd[ unique(readMatch[,1]) ] <- 1
    
##     # update results
    
##     fragSetCur <- vector( "list", length(existInd) )    
##     cur <- 1    
##     for ( j in 1:length(existInd) ) {
##         if ( existInd[j] == 1 ) {
##             fragSetCur[[cur]] <- 
##                 readCur[ matchList[[ as.character(j) ]], -1, drop=FALSE ]
##         } else {
##             fragSetCur[[cur]] <- matrix( NA )
##         }
##         cur <- cur + 1
##     }
    
##     # update name vector
    
##     nameVecCur <- 
##         paste( rep(chr,nrow(peakCur)), peakCur[,1], peakCur[,2], sep="_" )
    
##     # check proportion of forward reads for SET data
    
##     if ( PET == FALSE ) {
##         nFCur <- length(which( readCur[,4]=="F" ))
##         nAllCur <- nrow(readCur)
##     } else {
##         nFCur <- nAllCur <- NA
##     }
    
##     # calculate fragment length for PET data
    
##     if ( PET == TRUE ) {
##         fragLenCur <- readCur[,3] - readCur[,2] + 1
##     } else {
##         fragLenCur <- NA
##     }
    
##     return( list( fragSetCur=fragSetCur, nameVecCur=nameVecCur, 
##         peakChrCur=rep(chr,nrow(peakCur)),
##         peakStartCur=peakCur[,1], peakEndCur=peakCur[,2],
##         nFCur=nFCur, nAllCur=nAllCur, fragLenCur=fragLenCur ) )
## }

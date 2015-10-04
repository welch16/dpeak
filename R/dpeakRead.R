
# read bin-level data & process it

dpeakRead <- function( peakfile=NULL, readfile=NULL,
    fileFormat="eland_result", PET=FALSE, fragLen=200, 
    parallel=FALSE, nCore=8, tempDir=NULL, perl="perl" )
{            
    # process aligned read file
    
    if ( PET ) {
      message( "Info: Paired-end tag (PET) is assumed (PET=TRUE)." )
    } else {
      message( "Info: Single-end tag (SET) is assumed (PET=FALSE)." )
      message( "Info: Average fragment length is set as ",fragLen," (fragLen=",fragLen,")." )
    }

    #print( Sys.time() )
    
    message( "Info: Reading and processing aligned read file..." )
    
    if ( is.null(tempDir) ) {
      tempfileName <- tempfile( c("output","summary") )
    } else {
      tempfileName <- c( paste(tempDir,"output.txt",sep=""),
                        paste(tempDir,"summary.txt",sep="") )
    }

    # intermediate file name

    baseData <- .loadReadData(peakfile , readfile = readfile ,
        fileFormat = fileFormat,keepReads = TRUE,
        PET = PET,fragLen = fragLen, parallel = parallel,nCore = nCore, tempDir = tempDir,
        perl = perl)

                        
    ## .constructExtRead( infile=readfile, outfile=tempfileName[1],
    ##     summaryfile=tempfileName[2], fileFormat=fileFormat, 
    ##     PET=PET, fragLen=fragLen, perl=perl )            
    
    # read summary file (chrID, # lines)
    
    ## summaryInfo <- read.table( tempfileName[2], header=FALSE,
    ##     stringsAsFactors=FALSE, comment.char="", check.names=FALSE )
    ## colnames(summaryInfo) <- c("chrID","nline")   

    ## readChr <- summaryInfo[,1] 
    #print( Sys.time() )
        
    # process peak set & reads (assume BED file format)
    # - peak set: assume that first 3 columns of both files are chr, start, end
    
    ## message( "Info: Reading peak list..." )
    ## peakSet <- read.table( peakfile, header=FALSE, stringsAsFactors=FALSE )
    ## nPeak <- nrow(peakSet)
    ## gc()
    #print( Sys.time() )
    
    # match chromosomes between peak list and reads
    
    ##peakByChr <- split( peakSet[ , 2:3, drop=FALSE ], peakSet[,1] )
    ##peakChr <- names(peakByChr)
    
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


    fragSet <- baseData$fragSet
    numReads <- baseData$numReads
    stackedFragment <- baseData$stackedFragment
    seqDepth <- baseData$seqDepth

    idx <- which(numReads == 0)

    fragSet <- fragSet[-idx]
    numReads <- numReads[-idx]
    stackedFragment <- stackedFragment[-idx]
    
         
    ##fragSet <- vector( "list", nPeak )
    #peakChr <- peakStart <- peakEnd <- rep( NA, nPeak )


    nameVec <- names(fragSet)
    sep_names <- strsplit(nameVec,":")
    peakChr <- sapply(sep_names,function(x)x[1])
    sep_ranges <- strsplit(sapply(sep_names,function(x)x[2]),"-")

    peakStart <- sapply(sep_ranges,function(x)as.numeric(x[1]))
    peakEnd <- sapply(sep_ranges,function(x)as.numeric(x[2]))

    rm(sep_names,sep_ranges)

    nEmpty <- length(idx)
    if ( nEmpty == 0 ) {
      emptyList <- ""
    }else{
      emptyList <- names(idx)
    }
    
    
    ## nEmpty <- 0
    ## emptyList <- c()
    ## cur <- 1

    
    if ( PET == FALSE ) {
      nAll <- numReads
      if(parallel == TRUE){
        nF <- mclapply(fragSet,function(x)table(x$str)["F"],
            mc.cores = nCore)
      }else{
        nF <- lapply(fragSet,function(x)table(x$str)["F"])        
      }     
    }
    
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
    
    if ( PET == TRUE ) {    
        fragLen <- do.call(c,lapply(fragSet,width))
    }
        
    gc()    
    #print( Sys.time() )
        
    # check proportion of forward reads for SET data

    if ( PET == TRUE ) {
        Fratio <- 0.5
    } else {
      nF <- do.call(c,nF)
      Fratio <- sum(nF,na.rm = TRUE) / sum(nAll)
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
    
    ## nFrag <- unlist( lapply( fragSet, 
    ##     function(x) ifelse( !is.null(x), length(x), 0 )
    ##  ) )
    sumRead <- sum(numReads)
    medNumRead <- median(numReads)

    chrCommon <- unique(peakChr)
    nPeak <- length(fragSet)
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: Preprocessing summary\n" )
    cat( "------------------------------------------------------------\n" )
    cat( "Tag type: ",ifelse(PET,"PET","SET"),"\n", sep="" )
    cat( "Number of chromosomes: ",length(chrCommon),"\n", sep="" )
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

# match peak & reads

## .matchPeakRead <- function( chr, peakCur, outfileName, nRow, PET ) {
    
##     # read processed read file
##     # - PET read: chr, start, end
##     # - SET read: chr, position, strand, read length
    
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

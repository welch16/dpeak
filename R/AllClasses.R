
# Dpeak data class

setClass( Class="DpeakData",
    representation=representation(
        fragSet="list",
        
        PET="logical",
        fragLenTable="table",
        Fratio="numeric",
        aveFragLen="numeric",
        
        stackedFragment="list",
        
        peakChr="character",
        peakStart="numeric",
        peakEnd="numeric",
        
        emptyList="character"
    )
)

# Dpeak motif class

setClass( Class="DpeakMotif",
    representation=representation(
        motif="character",
        locMotif="list",
        
        peakChr="character",
        peakStart="numeric",
        peakEnd="numeric"
    )
)

# Dpeak fit class

setClass( Class="DpeakFit",
    representation=representation(
        fits="list",
        
        optFit="list",
        optMu="list",
        optPi="list",
        optPi0="list",
        optGamma="list",
        optDelta="list",
        optSigma="list",
        bicVec="list",
        aicVec="list",
        
        fragSet="list",
        PET="logical",
        fragLenTable="table",
        Fratio="numeric",
        aveFragLen="numeric",
        
        stackedFragment="list",
        
        peakChr="character",
        peakStart="numeric",
        peakEnd="numeric",
        
        estDeltaSigma="character",
        nTop="numeric",
        lbDelta="numeric",
        lbSigma="numeric",
        psize="numeric",
        maxComp="numeric", 
        pConst="numeric", 
        iterInit="numeric", 
        iterMain="numeric", 
        epsilon="numeric"
    )
)

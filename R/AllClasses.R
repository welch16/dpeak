
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
        
        estDelta="logical",
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

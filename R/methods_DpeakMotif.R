
# generic methods for "DpeakData" class

setMethod(
    f="show",
    signature="DpeakMotif",
    definition=function( object ) {    
      # extract objects
        
      motifIden <- object@motif
      locMotif <- object@locMotif
    	
      if ( all(is.na(unlist(locMotif))) ) {
      	# if no peak includes a motif
      	
      	cat( "------------------------------------------------------------\n" )
      	cat( "Info: Preprocessing summary\n" )
      	cat( "------------------------------------------------------------\n" )
      	cat( "Identified motif:\n", sep="" )
      	for ( i in 1:length(motifIden) ) {
      		cat( "\t",motifIden[i],"\n", sep="" )
      	}
      	cat( "Note: No peak contains detected motifs.\n")
      	cat( "------------------------------------------------------------\n" )
      } else {
		    # info about preprocessing
		    # - identified motif, # detected motifs
      	
      	listMotifVec <- list()
      	k <- 1
      	for ( i in 1:length(locMotif) ) {
      		if ( !is.na(locMotif[[i]][1]) ) {
      			listMotifVec[[k]] <- locMotif[[i]]
      			k <- k + 1
      		}
      	}
		    
		    cat( "------------------------------------------------------------\n" )
		    cat( "Info: Preprocessing summary\n" )
		    cat( "------------------------------------------------------------\n" )
		    cat( "Identified motif:\n", sep="" )
		    for ( i in 1:length(motifIden) ) {
		    	cat( "\t",motifIden[i],"\n", sep="" )
		    }
		    cat( "Number of peaks containing detected motifs: ",length(listMotifVec),"\n", sep="" )
		    cat( "Number of peaks missing detected motifs: ",( length(locMotif) - length(listMotifVec) ),"\n", sep="" )
		    cat( "Median number of regulatory sequences in each peak: ",median( sapply( listMotifVec, length ) ),"\n", sep="" )
		    cat( "\t(after excluding peaks missing detected motifs)\n" )
		    cat( "------------------------------------------------------------\n" )
    	}
    }
)

setMethod(
    f="print",
    signature="DpeakMotif",
    definition=function( x ) {
        warning( "'print' method for 'DpeakMotif' class is not supported yet." )
    }
)

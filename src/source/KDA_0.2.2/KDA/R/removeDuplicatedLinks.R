removeDuplicatedLinks <- function( linkpairs , directed = FALSE )
{
    if ( isTRUE( all.equal( dim( linkpairs )[1] , 1 ) ) )
	{
       return( linkpairs )
    }

    links <- paste( linkpairs[,1] , linkpairs[,2] , sep = "\t" )

    # 1. remove duplications 
    #
    cleanedlinkMatrix <- union( links , NULL )
# Why do we keep having these random printing things??
#	length( cleanedlinkMatrix )

    #ofname ="tmp.txt"
    #write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)

    # 2. remove inversed duplications
    #
    #linkMatrix <- read.delim(ofname, sep="\t", header=F) # not good for batch operation, ie, many of the same jobs running
    linkMatrix  <- getAllParts( cleanedlinkMatrix , "\t" )
# Again, another one!
#    dim(linkMatrix)
    linkMatrix <- as.matrix( linkMatrix )

    if ( directed )
	{
       return( linkMatrix )
    }

    if ( isTRUE( all.equal( dim( linkMatrix )[1] , 1 ) ) )
	{
		return( linkMatrix )
	}

    #  first, remove self-interactions are also removed
    #
    selSelfLinks <- linkMatrix[,1] == linkMatrix[,2]
    linkMatrix <- linkMatrix[!selSelfLinks,]
    cleanedlinkMatrix <- cleanedlinkMatrix[!selSelfLinks]

    reversedLinks <- paste( linkMatrix[,2] , linkMatrix[,1] , sep = "\t" )

    no.links <- length( reversedLinks )
    reversedLinksWithIndex <- cbind( reversedLinks , c( 1:no.links ) )

    merged <- merge( cleanedlinkMatrix , reversedLinksWithIndex , by.x = 1 , by.y = 1 , all = FALSE )
# What is it with these things?!?
#	dim(merged)

    if ( dim( merged )[1] > 0 )
	{
        merged <- as.matrix( merged )
        removedCols <- as.integer( merged[,2] )

        # cosntruct non-duplicated interactions
        #
        dupLinks <- cleanedlinkMatrix[removedCols]
        dupLinksRev <- reversedLinksWithIndex[removedCols]
        uniques <- NULL
        for ( i in c( 1:length( dupLinks ) ) )
        {
           found <- is.element( dupLinks[i] , uniques )
           foundRev <- is.element( dupLinksRev[i] , uniques )
           combined <- found | foundRev
           if ( !combined )
		   {
               uniques <- c( uniques , dupLinks[i] )
           }
        }
# ANOTHER ONE?!?
#		length(uniques)
        xlinkMatrix <- c( cleanedlinkMatrix[-removedCols] , uniques )
        #write.table( cleanedlinkMatrix[-removedCols], ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        #write.table( uniques, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE,append=T)
    }
	else
	{
        #write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        xlinkMatrix <- cleanedlinkMatrix
    }

    #linkMatrix <- read.delim(ofname, sep="\t", header=F)
    #dim(linkMatrix)
    #linkMatrix  = as.matrix(linkMatrix)

    return( getAllParts( xlinkMatrix , "\t" ) )
}


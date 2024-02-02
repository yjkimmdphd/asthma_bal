findNLayerNeighborsLinkPairs <- function( linkpairs , subnetNodes , nlayers = 1 , directed = FALSE )
{
  #linkpairs=linkpairs; subnetNodes=ineighbors; nlayers=nlayers-1; directed=directed
   merged <- merge( linkpairs , subnetNodes , by.x = 1 , by.y = 1 , all = FALSE )
   merged <- as.matrix( merged )

   if ( !directed )
   {
# undirected networks
     mergeleft <- merge( linkpairs , subnetNodes , by.x = 2 , by.y = 1 , all = FALSE )
     mergeleft <- as.matrix( mergeleft )
     mergeleft <- mergeleft[,c( 2 , 1)] # keep the original link direction
     merged <- rbind( merged , mergeleft )
   }

   if ( isTRUE( all.equal( dim( merged )[1] , 0 ) ) )
   {
	   return( NULL )
   }

   if ( isTRUE( all.equal( dim( merged )[1] , 1 ) ) &
	    isTRUE( all.equal( merged[1,1] , merged[1,2] ) ) )
   {
	   return( merged )
   }

   dim1 <- dim( merged )[1]
   if ( isTRUE( all.equal( dim1 , 0 ) ) )
   {
# no links
         return( NULL )
   }
   else if ( is.null( dim1 ) )
   {
# only one link
       merged <- rbind( merged )
   }
   merged <- removeDuplicatedLinks( merged , directed )
   if ( isTRUE( all.equal( nlayers , 1 ) ) )
   {
      return ( merged )
   }

   # get nodes   
   ineighbors <- union( merged[,1] , merged[,2] )

   if (nlayers==1){
      res=getSubnetwork_LinkPairs(linkpairs, subnetNodes=ineighbors)
      return (res)
   }


   # stop earlier if no change
   #
   common <- intersect( ineighbors , subnetNodes )
   if ( length( common ) == length( ineighbors ) )
   {
      return ( merged )
   }

   return ( findNLayerNeighborsLinkPairs( linkpairs , ineighbors , nlayers - 1 , directed ) )
}


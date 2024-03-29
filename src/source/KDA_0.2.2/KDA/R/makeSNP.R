makeSNP <- function( netpairsWtype , highlightNodes = NULL , edgecolorlevels , normColor = "grey" ,
		     highColor = "red" , normShape = "50" , highShape = "50", 
                     normNodeSize ="1",  highNodeSize="2", normFontSize ="1",    highFontSize="2",
                     directed = TRUE, legendtable = NA , snafile = "tmp.sna", kdaMatrix, featureData=NULL)
{
	if( ! is.null(snafile) ){ #added by MHW
    xfname = getFileName(snafile)
    fcys  = paste(xfname, "_cys.tsv", sep="")
    fcysn = paste(xfname, "_cys-nodes.tsv", sep="")
    write.table(netpairsWtype, fcys, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	} #added by MHW

    pcols = dim(netpairsWtype)[2]

    # consider both columns
    uniquenames    = union(as.character(netpairsWtype[,1]), as.character(netpairsWtype[,2]) )
    uniquenames    = sort(uniquenames)
    no.uniquenames = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    

    # 0. make link index: A1 A2 T & A I
    # A1 A2 T & A I ==> A1 A2 T I1
    #
    leftIdx = merge(netpairsWtype, name2idxMatrix, by.x=1, by.y=1, all=FALSE) 
    
    #  A1 A2 T I1 & A I ==> A2 A1 T I1 I2: I1 and I2 are the indices of A1 and A2 respectively
    #
    allIdx  = merge(leftIdx, name2idxMatrix, by.x=2, by.y=1, all=FALSE)

    no.pairs = dim(allIdx)[1]

    if(pcols==2){
      no.elevels = length(edgecolorlevels)
      greyIdx    = c(1:no.elevels) [edgecolorlevels=="grey"]
      linksIdxMatrix = cbind( allIdx[,c(3,4)], rep(greyIdx,no.pairs) )
    } else{
      linksIdxMatrix = allIdx[,c(4,5,3)]
    }


    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix 
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames),
                               rep(normNodeSize, no.uniquenames),
                               rep(normFontSize, no.uniquenames))
    } else{
      verticesMatrix = matrix("", no.uniquenames, 5)
      verticesMatrix[,1] = uniquenames
      verticesMatrix[,2] = rep(normColor,no.uniquenames)
      verticesMatrix[,3] = rep(normShape,no.uniquenames)
      verticesMatrix[,4] = rep(normNodeSize, no.uniquenames)
      verticesMatrix[,5] = rep(normFontSize, no.uniquenames)

      xcolor= rep(normColor,no.uniquenames);    names(xcolor) <- uniquenames
      xshape= rep(normShape,no.uniquenames);    names(xshape) <- uniquenames
      xNsize= rep(normNodeSize,no.uniquenames); names(xNsize) <- uniquenames
      xSsize= rep(normFontSize,no.uniquenames); names(xSsize) <- uniquenames

      # set color and shape for highlighted nodes
      #
      if ( !is.list(highlightNodes) ){
          highlightNodes2 = intersect(highlightNodes, uniquenames)
          xcolor[ highlightNodes2] = highColor[1]
          xshape[ highlightNodes2] = highShape[1]
          xNsize[ highlightNodes2] = highNodeSize[1]
          xSsize[ highlightNodes2] = highFontSize[1]
      }else{
         no.highnodeSets = length(highlightNodes)
         for(il in c(1:no.highnodeSets) ) {
             highlightNodes2 = intersect(highlightNodes[[il]], uniquenames)
             if(length(highlightNodes2)==0){next};
             xcolor[highlightNodes2] = highColor[il]
             xshape[highlightNodes2] = highShape[il]
             xNsize[highlightNodes2] = highNodeSize[il]
             xSsize[highlightNodes2] = highFontSize[il]
         }
      }
      verticesMatrix[,2] = as.character(xcolor)
      verticesMatrix[,3] = as.character(xshape)
      verticesMatrix[,4] = as.character(xNsize)
      verticesMatrix[,5] = as.character(xSsize)
    }

    #verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape", "size", "font_size")
	if(!is.null(featureData)){ #added by MHW
		verticesMatrix=data.frame(verticesMatrix,featureData[verticesMatrix[,'nodename'],],stringsAsFactors=FALSE)
		colnames(verticesMatrix)[5+1:ncol(featureData)]=colnames(featureData)
 	}

    #verticesMatrix2= verticesMatrix[,c(1:3)]
	verticesMatrix2= verticesMatrix[,-c(4,5)] #there could be more than 5 columns
    verticesMatrix2= merge(verticesMatrix2, kdaMatrix, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
    verticesMatrix3= as.matrix(verticesMatrix2)[,c(1:ncol(verticesMatrix2), ncol(verticesMatrix2)) ]
	colnames(verticesMatrix3)[ncol(verticesMatrix3)+c(-1,0)] <- c("size", "font_size") #edited by MHW
	if(is.null(snafile) ) return(list(verticesMatrix=verticesMatrix3)) #added by MHW
    write.table(verticesMatrix3, fcysn, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

    #**************************************************************************
    #
    # 3. output indexed netpairs
    #
    # Legend
 	snaConn=file(snafile,'w') #edited by MHW; replaced snafile by snaConn in the following
    if ( !is.null(legendtable) ) { #edited by MHW
      mhead = paste("Legend", dim(legendtable)[1], sep=" ")
      write.table(as.matrix(mhead), snaConn, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
      write.table(legendtable, snaConn, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
      # vertex
      write.table(as.matrix(header), snaConn, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
    } else{
      # vertex
      write.table(as.matrix(header), snaConn, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
    }

    write.table(verticesMatrix, snaConn, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

    #edge color
    write.table(t(as.matrix(c("edge_color_levels", edgecolorlevels) )), snaConn, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)

    # link pairs based index   
    #
    write.table(rbind(c("src","dst", "type")), snaConn, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
    write.table(linksIdxMatrix, snaConn, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
	close(snaConn)
    #return (c(fcys, fcysn))
	return (list(verticesMatrix=verticesMatrix,files=c(fcys, fcysn))) #edited by MHW
}

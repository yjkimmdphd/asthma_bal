\name{degree_ByLinkPairs}
\alias{degree_ByLinkPairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
degreeByLinkPairs(linkpairs, directed = F, cleangarbage = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linkpairs}{
%%     ~~Describe \code{linkpairs} here~~
}
  \item{directed}{
%%     ~~Describe \code{directed} here~~
}
  \item{cleangarbage}{
%%     ~~Describe \code{cleangarbage} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(linkpairs, directed=F, cleangarbage=F){

    codepair= c(0,1)  #[1] for no connection, [2] for connection

    edgesInNet = dim(linkpairs)[1]
    
    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(linkpairs[,1]) )
    allnodenames = c(allnodenames, as.character(linkpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)
    
    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    totallinks  = as.integer(nametable) # no of links for each node
    totalmatrix = cbind(names(nametable),  totallinks)

    if(directed){
        # outlines
        dnodenames = as.character(linkpairs[,1])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        outlinks = as.integer(as.matrix(iolinks[,2])) 


        # inlines
        dnodenames = as.character(linkpairs[,2])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        inlinks = as.integer(as.matrix(iolinks[,2])) 

    }else{
        inlinks  = totallinks
        outlinks = totallinks
    }

    #hubidx    = order(-totallinks)

    # output in/out links for each gene
    #
    linksMatrix            = cbind(inlinks, outlinks, totallinks)
    colnames(linksMatrix) <- c("inlinks", "outlinks", "totallinks")
    rownames(linksMatrix) <- uniquenames

    rm(inlinks)
    rm(outlinks)
    rm(totallinks)

    if(cleangarbage) {
        collect_garbage()
    }

    return ( data.frame(linksMatrix) )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

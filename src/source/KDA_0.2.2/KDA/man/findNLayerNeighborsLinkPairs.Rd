\name{findNLayerNeighbors_LinkPairs}
\alias{findNLayerNeighbors_LinkPairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findNLayerNeighbors_LinkPairs(linkpairs, subnetNodes, nlayers = 1, directed = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linkpairs}{
%%     ~~Describe \code{linkpairs} here~~
}
  \item{subnetNodes}{
%%     ~~Describe \code{subnetNodes} here~~
}
  \item{nlayers}{
%%     ~~Describe \code{nlayers} here~~
}
  \item{directed}{
%%     ~~Describe \code{directed} here~~
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
function(linkpairs, subnetNodes, nlayers=1, directed=F){

  #linkpairs=linkpairs; subnetNodes=ineighbors; nlayers=nlayers-1; directed=directed

   merged = merge(linkpairs, subnetNodes, by.x=1, by.y=1, all=F)
   merged <- as.matrix(merged)

   if ( dim(merged)[1] ==0){return (NULL) }

   if( (dim(merged)[1]==1) & (merged[1,1]==merged[1,2]) ){return(merged)}

   if(!directed) { # undirected networks
     mergeleft  = merge(linkpairs, subnetNodes, by.x=2, by.y=1, all=F)
     mergeleft  <- as.matrix(mergeleft)
     mergeleft  <- mergeleft[,c(2,1)] # keep the original link direction
     merged     <- rbind(merged, mergeleft)
   }

   dim1 = dim(merged)[1]
   if (dim1==0) { # no links
         return (NULL)
   }else if (is.null(dim1)){ # only one link
       merged = rbind(merged)
   }

   merged=removeDuplicatedLinks(merged, directed)
   
   if (nlayers==1){
      return (merged)
   }

   # get nodes   
   ineighbors =union(merged[,1], merged[,2])

   # stop earlier if no chnage
   #
   common = intersect(ineighbors, subnetNodes)
   if (length(common)==length(ineighbors)) {
      return (merged)
   }
   
   ret=findNLayerNeighbors_LinkPairs(linkpairs,ineighbors, nlayers-1, directed)

   return (ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

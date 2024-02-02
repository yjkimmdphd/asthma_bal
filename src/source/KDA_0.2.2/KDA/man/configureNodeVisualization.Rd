\name{configureNodeVisualization}
\alias{configureNodeVisualization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
configureNodeVisualization(allnodes, signature, kdaMatrix, bNodeSz=40, bFontSz=12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allnoes}{
%%     ~~Describe \code{allnoes} here~~
All the node names to be configured for visualization in cytoscape.
}
  \item{signature}{
%%     ~~Describe \code{signature} here~~
The input signature nodes/genes.
}
  \item{kdaMatrix}{
%%     ~~Describe \code{kdaMatrix} here~~
The result from key driver analysis
}

  \item{bNodeSz}{
%%     ~~Describe \code{bNodeSz} here~~
The default node size of the nodes
}

  \item{bFontSz}{
%%     ~~Describe \code{bFontSz} here~~
The default font size of the node labels}
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
Bin Zhang
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
function(allnodes, signature, kdaMatrix, bNodeSz=40, bFontSz=12) {

   # SIG--signature; NSIG--not signature; GKD--Global KeyDriver; LKD--Local KeyDriver; NKD--Not KeyDriver
   #
   xcategories = c("SIG_GKD", "SIG_LKD", "SIG_NKD", "NSIG_GKD", "NSIG_LKD", "NSIG_NKD"); xcat2=c("NSIG_GKD", "NSIG_LKD", "NSIG_NKD")
   xcolors     = c("red",     "blue",    "lightgreen",  "red",      "blue",     "grey");  names(xcolors)<- xcategories
   xshapes     = c("square",  "square",  "circle",   "circle",    "circle",   "circle");names(xshapes)<- xcategories
   xsizes      = c(3*bNodeSz, 2*bNodeSz, bNodeSz,   3*bNodeSz,  2*bNodeSz,  bNodeSz);     names(xsizes) <- xcategories
   xfontsz     = c(3*bFontSz, 2*bFontSz, bFontSz,   3*bFontSz,  2*bFontSz,  bFontSz);     names(xfontsz)<- xcategories

   no.nodes = length(allnodes)

   # legend table 
   legendtb = cbind(xcategories, xshapes, xcolors, xcolors, xsizes, xfontsz)
   colnames(legendtb) <- c("label", "shape", "color", "border", "node_size", "font_size")

   sigInNet = intersect(allnodes, signature)
   sig_status = rep("NSIG", no.nodes); names(sig_status) <- allnodes; sig_status[sigInNet]="SIG"
   kdr_status = rep("NKD",  no.nodes); names(kdr_status) <- allnodes; 

   nf.cols = dim(kdaMatrix)[2]; nf.rows = dim(kdaMatrix)[1]
   keydrvNames = NULL
   if(nf.rows>0) {
     keydrv  = as.integer(kdaMatrix[,nf.cols])

     # global driver
     keysel  = c(1:nf.rows)[keydrv==1]; 
     keydrvNames = kdaMatrix[keysel,1];
     kdr_status[keydrvNames] = "GKD"

     # local driver
     if(sum(keydrv==0)>0) {
        keysel  = c(1:nf.rows)[keydrv==0];
        keydrvNames = kdaMatrix[keysel,1];
        kdr_status[keydrvNames] = "LKD"
     }

     # combined signature-keydriver status
     #
     sigkdr_status=paste(sig_status, kdr_status, sep="_")
     hnList = tapply(allnodes, sigkdr_status, list) # make a list for each category
     sigkdr_names = names(hnList)

     isNonSig = intersect(xcat2, sigkdr_names) # if all nodes are signatures, we use only circle for display
     if(length(isNonSig)==0){
       xshapes     = c("circle",   "circle",   "circle",  "circle",   "circle",   "circle");names(xshapes)<- xcategories
     }

     # set up actual visualization properties
     yHighColor = xcolors[sigkdr_names]
     yHighShape = xshapes[sigkdr_names]
     yHighSize  = xsizes[sigkdr_names]
     yHighFontSZ= xfontsz[sigkdr_names]

   } else {
     hnList     = list(sigInNet) # highlight only signature
     yHighColor = c("brown")
     yHighShape = c("circle")
     yHighSize  = c("1")
     yHighFontSZ= c("1")
   }

   return( list(hnList, cbind(yHighColor, yHighShape, yHighSize, yHighFontSZ), legendtb) )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{keyDriverAnalysis}
\alias{keyDriverAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform key driver analysis on file.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Perform key driver analysis on an input tab-delimited network file with the first two columns representing the links.
}
\usage{
keyDriverAnalysis(inputnetwork, signature, directed=T, nlayer_expansion=1,
      nlayer_search=6, enrichedNodes_percent_cut=-1, boost_hubs=T, dynamic_search=T, 
      FET_pvalue_cut=0.05, use_corrected_pvalue=T, outputfile=NULL, expanded_network_as_signature=FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputnetwork}{
%%     ~~Describe \code{linkpairs} here~~
An input network file with the first 2-columns representing a links.
}
  \item{signature}{
%%     ~~Describe \code{signature} here~~
An input signature, i.e., a set of nodes.
}
  \item{directed}{
%%     ~~Describe \code{directed} here~~
A boolean variable indicating whether the input network is directed or not
}
  \item{nlayers}{
%%     ~~Describe \code{nlayers} here~~
The number of layers to be expanded
}
  \item{enrichedNodes_percent_cut}{
%%     ~~Describe \code{enrichedNodes_percent_cut} here~~
The minimal fraction of the nodes in a subnetwork overlapping the signature.
}
  \item{FET_pvalue_cut}{
%%     ~~Describe \code{FET_pvalue_cut} here~~
A threshold for Fisher's Exact Test pa value.
}
  \item{boost_hubs}{
%%     ~~Describe \code{boost_hubs} here~~
boost the nodes with high number of links as drivers or not.
}
  \item{dynamic_search}{
%%     ~~Describe \code{dynamic_search} here~~
Use dynamic search or not.
}
  \item{bonferroni_correction}{
%%     ~~Describe \code{bonferroni_correction} here~~
}

  \item{outputfile}{
%%     ~~Describe \code{outputfile} here~~
Output file name. By default the result is not output to files.
}
  \item{expanded_network_as_signature}{
%%     ~~Describe \code{expanded_network_as_signature} here~~
Use the expanded subnetwork as signature or not
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
\item{fkd}{Key driver status and statistis}
\item{cutmatrix}{Overall statistics}

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
function(inputnetwork, signature, directed=T, nlayer_expansion=1,
      nlayer_search=6, enrichedNodes_percent_cut=-1, boost_hubs=T, dynamic_search=T, 
      FET_pvalue_cut=0.05, use_corrected_pvalue=T, outputfile=NULL, expanded_network_as_signature=FALSE) 
{

   if (!is.null(outputfile)) {
      #onetFname   = paste(outputfile, outputDir, key2, ".pair", sep='')
      snpFname    = paste(outputfile, ".snp",  sep='')
      kdFname     = paste(outputfile, "_keydriver.xls",  sep='')
   }

   # overlap between network & signature
   wholenodes = union(inputnetwork[,1], inputnetwork[,2]); no.wholenodes=length(wholenodes)
   wholeOvlp  = intersect(wholenodes, signature); no.wholeOvlp = length(wholeOvlp)
   
   if(length(wholeOvlp)<=2) {return (NULL)}

   if(nlayer_expansion >=1 ) {
      # expand network by n-layer nearest neighbors
      expandNet = findNLayerNeighborsLinkPairs(linkpairs=inputnetwork, 
                  subnetNodes=signature, nlayers=nlayer_expansion, directed=directed)
   } else if(nlayer_expansion ==0 ){
      # no expansion
      expandNet = getSubnetworkLinkPairs(linkpairs=inputnetwork, subnetNodes=signature)
   } else{
      expandNet = inputnetwork
   }

   if(is.null(expandNet)) {return (NULL)}
   if(dim(expandNet)[1]<=10) {return (NULL)}

   dim(expandNet)
   print(paste("dim(expandNet): ", dim(expandNet)) )

   allnodes = sort(union(expandNet[,1], expandNet[,2])); no.nodes=length(allnodes)

   # convert IDs into indices
   netIdxSrc = getMatchedIndexFast(allnodes, expandNet[,1])
   netIdxDst = getMatchedIndexFast(allnodes, expandNet[,2])
   signatIdx = getMatchedIndexFast(allnodes, intersect(allnodes, signature))
   expandNetIdx = cbind(netIdxSrc, netIdxDst)

   ################################################################################################
   # 4. keydriver for a given network
   #
   #linkpairs=expandNetIdx; signature=signatIdx; background=c(no.wholenodes, no.wholeOvlp); directed=directed; nlayers=nlayer_search; enrichedNodes_percent_cut=enrichedNodes_percent_cut; FET_pvalue_cut=FET_pvalue_cut; boost_hubs=boost_hubs; dynamic_search=dynamic_search; bonferroni_correction=use_corrected_pvalue

   if (directed) {
     ret= keydriverInSubnetwork(linkpairs=expandNetIdx, signature=signatIdx, background=c(no.wholenodes, no.wholeOvlp),
              directed=directed, nlayers=nlayer_search, enrichedNodes_percent_cut=enrichedNodes_percent_cut, 
              FET_pvalue_cut=FET_pvalue_cut, 
              boost_hubs=boost_hubs, dynamic_search=dynamic_search, bonferroni_correction=use_corrected_pvalue, expanded_network_as_signature = expanded_network_as_signature)
   } else{
     ret= keydriverInSubnetwork(linkpairs=expandNetIdx, signature=signatIdx, background=c(no.wholenodes, no.wholeOvlp),
              directed=directed, nlayers=nlayer_search, enrichedNodes_percent_cut=enrichedNodes_percent_cut,
              FET_pvalue_cut=FET_pvalue_cut, 
              boost_hubs=boost_hubs, dynamic_search=dynamic_search, bonferroni_correction=use_corrected_pvalue, expanded_network_as_signature = expanded_network_as_signature)
   }

   if ( is.null(ret)) {return (NULL)}

   # retrieve results
   #
   fkd = ret[[1]]
   parameters = ret[[2]]

   fkd[,1] = allnodes[as.integer(fkd[,1])]   

   if (!is.null(outputfile)) {
      write.table(fkd, kdFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

      ################################################################################################
      #  output networks & key drivers for visualization
      #
      #     Cytoscape output: 1) network file - *_cys.txt 2) node property file: *_cys-nodes.txt
      #  * signature== genes need be corrected in other version
      # allnodes=allnodes; signature=signature; kdaMatrix=fkd; bNodeSz=40; bFontSz=12

      nodeprop = configureNodeVisualization(allnodes=allnodes, signature=signature, kdaMatrix=fkd)

      hnList     = nodeprop[[1]] # node subcategpries
      listprop   = nodeprop[[2]] # visual properties for each subcategory
      legend     = nodeprop[[3]] # legend table for visual propertie

      resf = makeSNP(netpairsWtype   = expandNet, 
               edgecolorlevels = c("grey"),
               highlightNodes  = hnList,
               normColor="grey",   highColor=listprop[,1],
               normShape="circle", highShape=listprop[,2],
               normNodeSize ="40",  highNodeSize =listprop[,3],
               normFontSize ="12",  highFontSize =listprop[,4],
               legendtable=legend, snafile=snpFname )

      result = list(expandNet, fkd, ret[[2]], getFileFullNameNopath(resf) )
      names(result) <- c("subnetwork", "keydrivers", "parameters", "files")
   } else{
      result = list(expandNet, fkd, ret[[2]])
      names(result) <- c("subnetwork", "keydrivers", "parameters")
   }

   return (result)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

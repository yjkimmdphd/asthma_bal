\name{getMatchedIndexFast}
\alias{getMatchedIndexFast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the matched index of a subvector
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMatchedIndexFast(cvector, subvect)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cvector}{
The first input vector
%%     ~~Describe \code{cvector} here~~
}
  \item{subvect}{
The second input vector that is supposed to be a  subset of the first one
%%     ~~Describe \code{subvect} here~~
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
function(cvector, subvect){
  fullindex = c(1:length(cvector) )
  orgIdx    = cbind(cvector, fullindex)

  index2    = c(1:length(subvect))
  subIdex   = cbind(subvect, index2)

  merged    = merge(subIdex, orgIdx, by.x=1, by.y=1, all.x=T)
  merged    = as.matrix(merged)

  if(dim(merged)[1]>1){
    od        = order(as.integer(merged[,2]))  # restore the original order of subvect
    merged    = merged[od, ]
  }
  
  outIndex  = as.integer(merged[,3])

  return (outIndex)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{concatenate}
\alias{concatenate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
concatenate(myvect, mysep = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myvect}{
%%     ~~Describe \code{myvect} here~~
}
  \item{mysep}{
%%     ~~Describe \code{mysep} here~~
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
function(myvect, mysep="")
{
  noitems = length(myvect)
  if (noitems==0){
    return ("")
  }else if (noitems==1){
    return (as.character(myvect) )
  }

  concatenated <- paste(as.character(myvect), sep="", collapse=mysep)
  return (concatenated)

  #tmpfn = "tmp.txt"
  #write.table(t(as.character(myvect)),tmpfn,sep=mysep,quote=FALSE, col.names=F, row.names=FALSE)
  #concatenated <- read.delim(tmpfn, sep="!", header=F)
  #return (as.character(as.matrix(concatenated) ))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

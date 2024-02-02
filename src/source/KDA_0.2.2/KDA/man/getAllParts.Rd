\name{getAllParts}
\alias{getAllParts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getAllParts(fullfnames, sep = "-", ret_len = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullfnames}{
%%     ~~Describe \code{fullfnames} here~~
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
}
  \item{ret_len}{
%%     ~~Describe \code{ret_len} here~~
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
function(fullfnames, sep="-", ret_len=F){

  splitted=unlist( strsplit(fullfnames[1], sep) )
  nn  = length(fullfnames)
  nf  = length(splitted)
  ret = matrix("", nn, nf)
  lens= rep(0, nn)
  for(i in c(1:nn) ) {
    each = fullfnames[i]
    splitted=unlist( strsplit(each, sep) )
    ino     = length(splitted)
    if(ino >=nf) {
       ret[i,] = splitted[1:nf]
    }else{
       ret[i,] = c(splitted, rep("",nf-ino ))
    }

    lens[i] = ino
  }

  if ( ret_len ){
     return (lens)
  } else {
     return (ret) 
  }
  
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

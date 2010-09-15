%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  cnr.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{cnr}
\alias{cnr}


\title{Simple creation of a CopyNumberRegions object}

\description{
 Simple creation of a CopyNumberRegions object containing a single region.
}

\usage{cnr(start, stop=NULL, mean=1, chromosome=1, xScale=1e+06, ...)}

\arguments{
  \item{start, stop}{Two \code{\link[base]{numeric}}s specifying the start and stop
    location of the region.  If \code{stop} is \code{\link[base]{NULL}}, then
    \code{start[2]} is used instead.}
  \item{mean}{A \code{\link[base]{numeric}} specifying the mean level of the region.}
  \item{chromosome}{An \code{\link[base]{integer}} specifying the chromosome ID.}
  \item{xScale}{The default scaling parameter used for \code{start}
    and \code{stop}, which is then also used for plotting genomic
    locations.}
  \item{...}{Additional arguments passed to
    \code{\link[aroma.core]{CopyNumberRegions}}.}
}

\value{
  Return an \code{\link[aroma.core]{CopyNumberRegions}} object.
}

\examples{
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build up a tumor CN profile
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pT <- cnr(1,1000, 2) +
     cnr(400,500) +
     cnr(600,800) +
     cnr(600,700) +
     cnr(100,200) - cnr(850,900);
print(pT);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate CN signals from this profile
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cn <- simulateRawCopyNumbers(pT, n=2000, sd=1/2);
print(cn);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, col="#aaaaaa", ylim=c(0,5));
drawLevels(pT, col="#ff0000", lwd=3);
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  callXXorXY.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{callXXorXY.numeric}
\alias{callXXorXY.numeric}
\alias{callXXorXY.numeric}

 \alias{callXXorXY}

 \title{Calls XX or XY from ChrX allele B fractions of a normal sample}

 \description{
   Calls XX or XY from ChrX allele B fractions of a normal sample.
 }

 \usage{\method{callXXorXY}{numeric}(betaX, betaY=NULL, flavor=c("density"), adjust=1.5, ..., censorAt=c(-0.5, +1.5), verbose=FALSE)}

 \arguments{
  \item{betaX}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} containing ChrX allele B fractions.}
  \item{betaY}{A optional \code{\link[base]{numeric}} \code{\link[base]{vector}} containing ChrY allele B fractions.}
  \item{flavor}{A \code{\link[base]{character}} string specifying the type of algorithm used.}
  \item{adjust}{A postive \code{\link[base]{double}} specifying the amount smoothing for
    the empirical density estimator.}
  \item{...}{Additional arguments passed to
    \code{\link[aroma.light]{findPeaksAndValleys}}.}
  \item{censorAt}{A \code{\link[base]{double}} \code{\link[base]{vector}} of length two specifying the range
    for which values are considered finite.  Values below (above) this
    range are treated as -\code{\link[base:is.finite]{Inf}} (+\code{\link[base:is.finite]{Inf}}).}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
 }

 \value{
   Returns a ...
 }

 \section{Missing and non-finite values}{
   Missing and non-finite values are dropped before trying to call XX or XY.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Internally \code{\link[aroma.light]{findPeaksAndValleys}} is used to identify
   the thresholds.
 }
\keyword{methods}

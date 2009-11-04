%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  callNaiveGenotypes.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{callNaiveGenotypes.numeric}
\alias{callNaiveGenotypes.numeric}
\alias{callNaiveGenotypes.numeric}

 \alias{callNaiveGenotypes}

 \title{Calls genotypes in a normal sample}

 \description{
   Calls genotypes in a normal sample.
 }

 \usage{\method{callNaiveGenotypes}{numeric}(y, cn=rep(2, length(y)), flavor=c("density"), adjust=1.5, ..., censorAt=c(-0.5, +1.5), verbose=FALSE)}

 \arguments{
  \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J containing allele B fractions
    for a normal sample.}
  \item{cn}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J specifying the true
    total copy number in \eqn{\{0,1,2,NA\}} at each locus.  This can be
    used to specify which loci are diploid and which are not, e.g.
    autosomal and sex chromosome copy numbers.}
  \item{flavor}{A \code{\link[base]{character}} string specifying the type of algorithm used.}
  \item{adjust}{A postive \code{\link[base]{double}} specifying the amount smoothing for
    the empirical density estimator.}
  \item{...}{Additional arguments passed to \code{\link{findPeaksAndValleys}}().}
  \item{censorAt}{A \code{\link[base]{double}} \code{\link[base]{vector}} of length two specifying the range
    for which values are considered finite.  Values below (above) this
    range are treated as -\code{\link[base:is.finite]{Inf}} (+@inf).}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
 }

 \value{
   Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J containing the genotype calls
   in allele B fraction space, that is, in [0,1] where 1/2 corresponds
   to a heterozygous call, and 0 and 1 corresponds to homozygous A
   and B, respectively.
   Non called genotypes have value \code{\link[base]{NA}}.
 }

 \examples{
subplots(3, ncol=1)
par(mar=c(2,4,4,1)+0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A bimodal distribution
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAA <- rnorm(n=10000, mean=0, sd=0.1)
xBB <- rnorm(n=10000, mean=1, sd=0.1)
x <- c(xAA,xBB)
fit <- findPeaksAndValleys(x)
print(fit)
calls <- callNaiveGenotypes(x, cn=rep(1,length(x)), verbose=-20)
xc <- split(x, calls)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="(AA,BB)")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution with missing values
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAB <- rnorm(n=10000, mean=1/2, sd=0.1)
x <- c(xAA,xAB,xBB)
x[sample(length(x), size=0.05*length(x))] <- NA;
x[sample(length(x), size=0.01*length(x))] <- -Inf;
x[sample(length(x), size=0.01*length(x))] <- +Inf;
fit <- findPeaksAndValleys(x)
print(fit)
calls <- callNaiveGenotypes(x)
xc <- split(x, calls)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="(AA,AB,BB)")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution with clear separation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAA <- rnorm(n=10000, mean=0, sd=0.02)
xAB <- rnorm(n=10000, mean=1/2, sd=0.02)
xBB <- rnorm(n=10000, mean=1, sd=0.02)
x <- c(xAA,xAB,xBB)
fit <- findPeaksAndValleys(x)
print(fit)
calls <- callNaiveGenotypes(x)
xc <- split(x, calls)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="(AA',AB',BB')")
abline(v=fit$x)
}

 \section{Missing and non-finite values}{
   A missing value always gives a missing (\code{\link[base]{NA}}) genotype call.
   Negative infinity (-\code{\link[base:is.finite]{Inf}}) always gives genotype call 0.
   Positive infinity (+\code{\link[base:is.finite]{Inf}}) alwsys gives genotype call 1.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Internally \code{\link{findPeaksAndValleys}}() is used to identify the thresholds.
 }
\keyword{methods}

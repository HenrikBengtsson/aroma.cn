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

 \usage{\method{callNaiveGenotypes}{numeric}(y, cn=rep(2, length(y)), flavor=c("density"), adjust=1.5, ...)}

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
x1 <- rnorm(n=10000, mean=0, sd=1)
x2 <- rnorm(n=10000, mean=4, sd=1)
x <- c(x1,x2)
fit <- findPeaksAndValleys(x)
calls <- callNaiveGenotypes(x, cn=rep(1,length(x)))
xc <- split(x, calls)
print(fit)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="c(x1,x2)")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x3 <- rnorm(n=10000, mean=8, sd=1)
x <- c(x1,x2,x3)
fit <- findPeaksAndValleys(x)
calls <- callNaiveGenotypes(x)
xc <- split(x, calls)
print(fit)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="c(x1,x2,x3)")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution with clear separation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x1b <- rnorm(n=10000, mean=0, sd=0.1)
x2b <- rnorm(n=10000, mean=4, sd=0.1)
x3b <- rnorm(n=10000, mean=8, sd=0.1)
x <- c(x1b,x2b,x3b)
fit <- findPeaksAndValleys(x, tol=0)
calls <- callNaiveGenotypes(x, tol=0)
xc <- split(x, calls)
print(fit)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="c(x1b,x2b,x3b)")
abline(v=fit$x)
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Internally \code{\link{findPeaksAndValleys}}() is used to identify the thresholds.
 }
\keyword{methods}

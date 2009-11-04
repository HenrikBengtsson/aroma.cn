%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findPeaksAndValleys.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{findPeaksAndValleys.numeric}
\alias{findPeaksAndValleys.numeric}
\alias{findPeaksAndValleys.numeric}

 \alias{findPeaksAndValleys}

 \title{Finds extreme points in the empirical density estimated from data}

 \description{
   Finds extreme points in the empirical density estimated from data.
 }

 \usage{\method{findPeaksAndValleys}{numeric}(x, na.rm=TRUE, ..., tol=0)}

 \arguments{
  \item{x}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} containing data points.}
  \item{na.rm}{If \code{\link[base:logical]{TRUE}}, missing values are dropped, otherwise not.}
  \item{...}{Arguments passed to \code{\link[stats]{density}}.}
  \item{tol}{A non-negative \code{\link[base]{numeric}} threshold specifying the minimum
    density at the extreme point in order to accept it.}
 }

 \value{
   Returns a @dataframe of "peaks" and "valleys" filtered by \code{tol}.
 }

 \examples{
subplots(3, ncol=1)
par(mar=c(2,4,4,1)+0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A unimodal distribution
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x1 <- rnorm(n=10000, mean=0, sd=1)
x <- x1
fit <- findPeaksAndValleys(x)
print(fit)
plot(density(x), lwd=2, main="x1")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x2 <- rnorm(n=10000, mean=4, sd=1)
x3 <- rnorm(n=10000, mean=8, sd=1)
x <- c(x1,x2,x3)
fit <- findPeaksAndValleys(x)
print(fit)
plot(density(x), lwd=2, main="c(x1,x2,x3)")
abline(v=fit$x)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution with clear separation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x1b <- rnorm(n=10000, mean=0, sd=0.1)
x2b <- rnorm(n=10000, mean=4, sd=0.1)
x3b <- rnorm(n=10000, mean=8, sd=0.1)
x <- c(x1b,x2b,x3b)
fit <- findPeaksAndValleys(x, tol=0)
print(fit)
plot(density(x), lwd=2, main="c(x1b,x2b,x3b)")
abline(v=fit$x)
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{methods}

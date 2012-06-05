%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizePrincipalCurve.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizePrincipalCurve.matrix}
\alias{normalizePrincipalCurve.matrix}
\alias{normalizePrincipalCurve.matrix}


\encoding{latin1}

\title{Fit a principal curve in K dimensions}

\description{
  Fit a principal curve in K dimensions.
}

\usage{\method{normalizePrincipalCurve}{matrix}(x, ..., center=TRUE, returnFit=FALSE)}

\arguments{
 \item{x}{An NxK \code{\link[base]{matrix}} (K>=2) where the columns represent the dimension.}
 \item{...}{Additional arguments passed to
     \code{\link[aroma.light]{fitPrincipalCurve}} used for fitting the model.}
 \item{center}{If \code{\link[base:logical]{TRUE}}, normalized data is centered such that the median
     signal in each dimension is at zero.}
 \item{returnFit}{If \code{\link[base:logical]{TRUE}}, the fitted principal curve parameters are
     returned as an attribute.}
}

\value{
  Returns an NxK \code{\link[base]{matrix}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] Hastie, T. and Stuetzle, W, \emph{Principal Curves}, JASA, 1989.
}

\seealso{
  \code{\link[aroma.light]{fitPrincipalCurve}} and
  \code{\link[aroma.light]{backtransformPrincipalCurve}}.
}
\keyword{methods}

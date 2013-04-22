%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizePrincipalCurve.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizePrincipalCurve.matrix}
\alias{normalizePrincipalCurve.matrix}

\encoding{latin1}

\title{Normalizes data in K dimensions using principal curves}

\description{
  Normalizes data in K dimensions using principal curves such that afterward the data cluster (approximately
  linearly) along the diagonal (in K dimensions).
}

\usage{
\method{normalizePrincipalCurve}{matrix}(x, ..., center=TRUE, returnFit=FALSE)
}

\arguments{
 \item{x}{An NxK \code{\link[base]{matrix}} where the columns represent the
     (K >= 2) dimensions.}
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

\references{
  [1] Hastie, T. and Stuetzle, W, \emph{Principal Curves}, JASA, 1989.
}

\seealso{
  \code{\link[aroma.light]{fitPrincipalCurve}} and
  \code{\link[aroma.light]{backtransformPrincipalCurve}}.
}

\author{Henrik Bengtsson}


\keyword{methods}
\keyword{internal}

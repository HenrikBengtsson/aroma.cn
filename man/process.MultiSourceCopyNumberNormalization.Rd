%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiSourceCopyNumberNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{process.MultiSourceCopyNumberNormalization}
\alias{process.MultiSourceCopyNumberNormalization}
\alias{MultiSourceCopyNumberNormalization.process}
\alias{process.MultiSourceCopyNumberNormalization}
\alias{process,MultiSourceCopyNumberNormalization-method}


\title{Normalizes all samples}

\description{
 Normalizes all samples.
}

\usage{\method{process}{MultiSourceCopyNumberNormalization}(this, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{list}} of K \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}}:s.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.cn:fit.MultiSourceCopyNumberNormalization]{*fit}()}.
  For more information see \code{\link{MultiSourceCopyNumberNormalization}}.
}
\keyword{internal}
\keyword{methods}

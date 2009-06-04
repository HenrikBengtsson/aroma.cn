%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiSourceCopyNumberNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fitOne.MultiSourceCopyNumberNormalization}
\alias{fitOne.MultiSourceCopyNumberNormalization}
\alias{MultiSourceCopyNumberNormalization.fitOne}
\alias{fitOne.MultiSourceCopyNumberNormalization}
\alias{fitOne,MultiSourceCopyNumberNormalization-method}


\title{Fits the multi-source model for one sample}

\description{
 Fits the multi-source model for one sample.
 The model is fitted on the subset of units returned
 by \code{\link[aroma.cn:getSubsetToFit.MultiSourceCopyNumberNormalization]{*getSubsetToFit}()}.
}

\usage{\method{fitOne}{MultiSourceCopyNumberNormalization}(this, dfList, ..., force=FALSE, .retData=FALSE, verbose=FALSE)}

\arguments{
  \item{name}{A \code{\link[base]{character}} string specifying the sample name of interest.}
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{FALSE}}, cached model fits are returned, otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{list}} of transforms.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  This method is called internally by \code{\link[aroma.cn:fit.MultiSourceCopyNumberNormalization]{*fit}()}.
  For more information see \code{\link{MultiSourceCopyNumberNormalization}}.
}
\keyword{internal}
\keyword{methods}

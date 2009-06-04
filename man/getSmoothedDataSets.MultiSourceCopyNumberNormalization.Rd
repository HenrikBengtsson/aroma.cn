%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiSourceCopyNumberNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSmoothedDataSets.MultiSourceCopyNumberNormalization}
\alias{getSmoothedDataSets.MultiSourceCopyNumberNormalization}
\alias{MultiSourceCopyNumberNormalization.getSmoothedDataSets}
\alias{getSmoothedDataSets.MultiSourceCopyNumberNormalization}
\alias{getSmoothedDataSets,MultiSourceCopyNumberNormalization-method}


\title{Gets the data sets smoothed toward the UGP file}

\description{
 Gets the data sets smoothed toward the UGP file.
}

\usage{\method{getSmoothedDataSets}{MultiSourceCopyNumberNormalization}(this, ..., verbose=FALSE)}

\arguments{
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{list}} of K \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}}:s.
}

\details{
  This method smooth each data set, each array, and each chromosome
  toward the target (smoothing) UGP file independently of everything else.

  The resulting data sets are stored in a separate location where they
  will be located automatically in future sessions.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.cn:getFitAromaUgpFile.MultiSourceCopyNumberNormalization]{*getFitAromaUgpFile}()}.
  For more information see \code{\link{MultiSourceCopyNumberNormalization}}.
}
\keyword{internal}
\keyword{methods}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiSourceCopyNumberNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{extractTupleOfDataFiles.MultiSourceCopyNumberNormalization}
\alias{extractTupleOfDataFiles.MultiSourceCopyNumberNormalization}
\alias{MultiSourceCopyNumberNormalization.extractTupleOfDataFiles}
\alias{extractTupleOfDataFiles.MultiSourceCopyNumberNormalization}
\alias{extractTupleOfDataFiles,MultiSourceCopyNumberNormalization-method}


\title{Gets a list of data files for a particular name across several data sets}

\description{
 Gets a list of data files for a particular name across several data sets.
}

\usage{\method{extractTupleOfDataFiles}{MultiSourceCopyNumberNormalization}(this, dsList, name, ..., na.rm=FALSE, verbose=FALSE)}

\arguments{
  \item{name}{A \code{\link[base]{character}} string specifying the sample name of interest.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{list}} of K \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}}:s.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MultiSourceCopyNumberNormalization}}.
}
\keyword{internal}
\keyword{methods}

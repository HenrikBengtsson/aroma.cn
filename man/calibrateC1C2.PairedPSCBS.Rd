%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.calibrateC1C2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{calibrateC1C2.PairedPSCBS}
\alias{calibrateC1C2.PairedPSCBS}
\alias{PairedPSCBS.calibrateC1C2}
\alias{calibrateC1C2.PairedPSCBS}
\alias{calibrateC1C2,PairedPSCBS-method}


\title{Calibrates ASCN signals in (C1,C2) space based on region-based PSCN estimates}

\description{
 Calibrates ASCN signals in (C1,C2) space based on region-based PSCN estimates as given by the PSCBS segmentation method.
}

\usage{\method{calibrateC1C2}{PairedPSCBS}(fit, ..., force=FALSE, cache=TRUE, verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[psCBS]{segmentByPairedPSCBS}}.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a calibrated PairedPSCBS fit object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}

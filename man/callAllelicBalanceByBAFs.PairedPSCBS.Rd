%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.CALL2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callAllelicBalanceByBAFs.PairedPSCBS}
\alias{callAllelicBalanceByBAFs.PairedPSCBS}
\alias{PairedPSCBS.callAllelicBalanceByBAFs}
\alias{callAllelicBalanceByBAFs.PairedPSCBS}
\alias{callAllelicBalanceByBAFs,PairedPSCBS-method}


\title{Calls regions that are in allelic balance}

\description{
 Calls regions that are in allelic balance from the allele B fractions (BAF).
}

\usage{\method{callAllelicBalanceByBAFs}{PairedPSCBS}(fit, maxScore="auto", ..., force=FALSE, cache=FALSE, verbose=FALSE)}

\arguments{
  \item{fit}{A \code{\link[PSCBS]{PairedPSCBS}} fit object as returned by
    \code{\link[PSCBS]{segmentByPairedPSCBS}}.}
  \item{maxScore}{A positive \code{\link[base]{double}} threshold.
    If \code{"auto"}, the threshold is estimated empirically.}
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, an already called object is skipped, otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link[PSCBS]{PairedPSCBS}} fit object
  where columns for allelic imbalance scores and
  p-values as well as allelic balance calls are added.
}

\examples{

if (require("PSCBS") && Sys.getenv("_R_CHECK_FULL_") != "") {

# Load example ASCN data
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS");
data <- R.utils::loadObject(pathname);

# PSCBS segmentation
fit <- segmentByPairedPSCBS(data, verbose=-10);

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

fitN2 <- callAllelicBalanceByBAFs(fitN, verbose=-10);
print(fitN2);

} # if (require("PSCBS"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, \code{\link{testAllelicBalanceByBAFs}}() is used.

  Note that this AB caller differs from the default one in the
  \pkg{PSCBS} package, cf. \code{\link[PSCBS]{callAB.PairedPSCBS}}.
}


\keyword{internal}
\keyword{methods}

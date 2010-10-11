%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.EXTS2.R
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

\usage{\method{callAllelicBalanceByBAFs}{PairedPSCBS}(fit, maxScore=4, ..., force=FALSE, cache=TRUE, verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[psCBS]{segmentByPairedPSCBS}}.}
  \item{maxScore}{A positive \code{\link[base]{double}} threshold.}
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, an already called object is skipped, otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a PairedPSCBS fit object where columns for
  allelic imbalance scores and p-values as well as allelic
  balance calls are added.
}

\examples{
if (require("psCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
data <- R.utils::loadObject(pathname);
R.oo::attachLocally(data);
x <- position;

# PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, x=x, verbose=-10);

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

fitN2 <- callAllelicBalanceByBAFs(fitN, verbose=-10);
print(fitN2);

} # if (require("psCBS"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}

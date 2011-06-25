%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.EXTS2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{callCopyNeutralRegions.PairedPSCBS}
\alias{callCopyNeutralRegions.PairedPSCBS}
\alias{PairedPSCBS.callCopyNeutralRegions}
\alias{callCopyNeutralRegions.PairedPSCBS}
\alias{callCopyNeutralRegions,PairedPSCBS-method}


\title{Calls regions that are copy neutral}

\description{
 Calls regions that are copy neutral from the allele B fractions (BAF).
}

\usage{\method{callCopyNeutralRegions}{PairedPSCBS}(fit, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[PSCBS]{segmentByPairedPSCBS}}.}
  \item{...}{Additional arguments passed to
    \code{\link[aroma.cn]{findNeutralCopyNumberState}}.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, an already called object is skipped, otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a PairedPSCBS fit object where a column with the copy-neutral call.
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

fitN2 <- callCopyNeutralRegions(fitN, verbose=-10);
print(fitN2);

} # if (require("psCBS"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}

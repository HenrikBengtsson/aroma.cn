%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.CALL2.R
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

\usage{\method{callCopyNeutralRegions}{PairedPSCBS}(fit, ...)}

\arguments{
  \item{fit}{A \code{\link[PSCBS]{PairedPSCBS}} fit object as returned by
    \code{\link[PSCBS]{segmentByPairedPSCBS}}.}
  \item{...}{Additional arguments passed to
    \code{\link[aroma.cn]{findNeutralCopyNumberState}} and
    \code{\link[PSCBS]{callCopyNeutral.PairedPSCBS}}.}
}

\value{
  Returns a \code{\link[PSCBS]{PairedPSCBS}} fit object
  where a column with the copy-neutral call.
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

fitN2 <- callCopyNeutralRegions(fitN, verbose=-10);
print(fitN2);

} # if (require("PSCBS"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}

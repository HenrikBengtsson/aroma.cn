%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeBAFsByRegions.PairedPSCBS}
\alias{normalizeBAFsByRegions.PairedPSCBS}
\alias{PairedPSCBS.normalizeBAFsByRegions}
\alias{normalizeBAFsByRegions.PairedPSCBS}
\alias{normalizeBAFsByRegions,PairedPSCBS-method}


\title{Normalizes allele B fractions (BAFs) based on region-based PSCN estimates}

\description{
 Normalizes allele B fractions (BAFs) based on region-based PSCN estimates as given by the PSCBS segmentation method.
}

\usage{\method{normalizeBAFsByRegions}{PairedPSCBS}(fit, by=c("betaTN", "betaT"), ..., force=FALSE, cache=TRUE, verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[PSCBS]{segmentByPairedPSCBS}}.}
  \item{by}{A \code{\link[base]{character}} string specifying if the normalization function
    should be estimated based on TumorBoost normalized or non-normalized
    tumor allele B fractions (BAFs).}
  \item{...}{Additional arguments passed
    \code{\link[aroma.cn]{normalizeMirroredBAFsByRegions}}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a PairedPSCBS fit object where the region-level
  decrease-in-heterozygosity (DH) means have been normalized,
  as well as the locus-specific tumor allele B fractions.
}

\details{
  Note that his normalization method depends on the segmentation
  results. Hence, it recommended \emph{not} to resegment the
  normalized signals returned by this, because such a segmentation
  will be highly dependent on the initial segmentation round.
}

\examples{
if (require("PSCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS");
data <- R.utils::loadObject(pathname);
R.oo::attachLocally(data);
x <- position;

# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

# PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, x=x, verbose=-10);

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

devSet("tracks");
subplots(10, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plot(fit, subplots=TRUE);
plot(fitN, subplots=TRUE);


devSet("C1C2");
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);
plotC1C2(fit, Clim=Clim);
linesC1C2(fit);
title(main="(C1,C2)");

plotC1C2(fitN, Clim=Clim);
linesC1C2(fitN);
title(main="(C1,C2) - adjusted");

plotC1C2(fit, col="gray", Clim=Clim);
linesC1C2(fit, col="gray");
pointsC1C2(fitN);
linesC1C2(fitN);
title(main="(C1,C2) - both");

} # if (require("PSCBS"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[aroma.cn]{normalizeMirroredBAFsByRegions}} is used.
}


\keyword{internal}
\keyword{methods}

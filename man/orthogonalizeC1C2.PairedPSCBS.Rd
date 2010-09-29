%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{orthogonalizeC1C2.PairedPSCBS}
\alias{orthogonalizeC1C2.PairedPSCBS}
\alias{PairedPSCBS.orthogonalizeC1C2}
\alias{orthogonalizeC1C2.PairedPSCBS}
\alias{orthogonalizeC1C2,PairedPSCBS-method}


\title{Orthogonalizes (C1,C2) based on region-based PSCN estimates}

\description{
 Orthogonalizes (C1,C2) based on region-based PSCN estimates as given by the PSCBS segmentation method.
}

\usage{\method{orthogonalizeC1C2}{PairedPSCBS}(fit, ..., debugPlot=TRUE, verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[psCBS]{segmentByPairedPSCBS}}.}
  \item{...}{Additional arguments passed
    \code{\link[aroma.light]{findPeaksAndValleys}}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a PairedPSCBS fit object.
}

\examples{
if (require("psCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
data <- R.utils::loadObject(pathname);
R.oo::attachLocally(data);
x <- position;

# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

# PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, x=x, verbose=-10);
fit <- postsegmentTCN(fit);
print(fit);

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

# Orthogonalize
fitO <- orthogonalizeC1C2(fitN);


devSet("tracks");
subplots(10, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plot(fitN, add=TRUE);
plot(fitO, add=TRUE);

devSet("C1C2");
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);
plotC1C2(fitN, Clim=Clim);
linesC1C2(fitN);
title(main="(C1,C2) - adjusted");

plotC1C2(fitO, Clim=Clim);
linesC1C2(fitO);
title(main="(C1,C2) - orthogonalized");

plotC1C2(fitN, col="gray", Clim=Clim);
linesC1C2(fitN, col="gray");
pointsC1C2(fitO);
linesC1C2(fitO);
title(main="(C1,C2) - both");

} # if (require("psCBS"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}

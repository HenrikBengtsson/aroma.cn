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

\usage{\method{calibrateC1C2}{PairedPSCBS}(fit, ..., force=FALSE, cache=TRUE, debug=FALSE, verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[psCBS]{segmentByPairedPSCBS}}.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a calibrated PairedPSCBS fit object.
}

\examples{\dontrun{
 if (require("psCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
# pathname <- "wholeGenomeData/GSE20584,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY/GenomeWideSNP_6/GSE20584,GSM517071vGSM517072,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY.xdr";
data <- R.utils::loadObject(pathname);
R.oo::attachLocally(data);
x <- position;

# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

# PSCBS segmentation
cc <- 1;
ww <- which(data$chromosome == cc);

fit <- segmentByPairedPSCBS(CT[ww], betaT=betaT[ww], betaN=betaN[ww], x=x[ww], chromosome=cc, verbose=-10);
fit <- postsegmentTCN(fit);

# Calibrate in (C1,C2) space
fitC <- calibrateC1C2(fit, verbose=verbose);

devSet("tracks");
subplots(10, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plot(fit, add=TRUE);
plot(fitC, add=TRUE);

devSet("C1C2");
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);
plotC1C2(fit, Clim=Clim);
linesC1C2(fit);
title(main="(C1,C2) - raw");

plotC1C2(fitC, Clim=Clim);
linesC1C2(fitC);
title(main="(C1,C2) - calibrated");

} # if (require("psCBS"))

}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}

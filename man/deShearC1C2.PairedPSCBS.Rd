%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.deShearC1C2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{deShearC1C2.PairedPSCBS}
\alias{deShearC1C2.PairedPSCBS}
\alias{PairedPSCBS.deShearC1C2}
\alias{deShearC1C2.PairedPSCBS}
\alias{deShearC1C2,PairedPSCBS-method}


\title{Correct for shearing in (C1,C2) space based on region-based PSCN estimates}

\description{
 Correct for shearing in (C1,C2) space based on region-based PSCN estimates as given by the PSCBS segmentation method.
}

\usage{\method{deShearC1C2}{PairedPSCBS}(fit, ..., dirs=c("|-", "-", "|", "X", "|,-", "-,|", "|-,X", "|,-,X", "-,|,X"), verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[PSCBS]{segmentByPairedPSCBS}}.}
  \item{adjust}{A \code{\link[base]{numeric}} adjusting the bandwidth of the empirical
    density estimator of line (changepoint) directions.}
  \item{weightFlavor}{A \code{\link[base]{character}} string specifying how weights are
    generated for lines (changepoints).}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a PairedPSCBS fit object.
}

\examples{

if (require("PSCBS") && Sys.getenv("_R_CHECK_FULL_") != "") {

# Load example ASCN data
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS");
# pathname <- "wholeGenomeData/GSE20584,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY/GenomeWideSNP_6/GSE20584,GSM517071vGSM517072,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY.xdr";
data <- R.utils::loadObject(pathname);
R.oo::attachLocally(data);

# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

# PSCBS segmentation
cc <- 1;
ww <- which(data$chromosome == cc);

fit <- segmentByPairedPSCBS(CT[ww], betaT=betaT[ww], betaN=betaN[ww], x=x[ww], chromosome=cc, verbose=-10);
fit <- postsegmentTCN(fit);

## Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

# Call allelic balance
fitC <- callAllelicBalanceByBAFs(fitN);
ww <- which(fitC$output$abCall);
fitC$output[ww, "dhMean"] <- 0;

# deShear
fitD <- deShearC1C2(fitC, verbose=verbose);

devSet("tracks");
subplots(2*3, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plotTracks(fit, subplots=FALSE);
plotTracks(fitD, subplots=FALSE);

devSet("C1C2");
Clim <- c(0,4);
plotC1C2(fitC, Clim=Clim);
linesC1C2(fitC);
title(main="(C1,C2) - adjusted");

devSet("C1C2,deShear");
Clim <- c(0,4);
plotC1C2(fitD, Clim=Clim);
linesC1C2(fitD);
title(main="(C1,C2) - deSheared");

} # if (require("PSCBS"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}

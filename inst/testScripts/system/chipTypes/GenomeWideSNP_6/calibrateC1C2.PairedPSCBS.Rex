library("aroma.cn");
library("PSCBS");
#sourceTo("../aroma.cn/R/PairedPSCBS.NORM.R", modifiedOnly=TRUE);


verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load ASCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "wholeGenomeData";

dataSet <- "GSE20584,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";
chipType <- "GenomeWideSNP_6";
filename <- "GSE20584,GSM517071vGSM517072,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY.xdr";

path <- file.path(rootPath, dataSet, chipType);
path <- ".";
pathname <- Arguments$getReadablePathname(filename, path=path);

if (!exists("data", mode="list")) {
  verbose && enter(verbose, "Loading ASCN data");
  dataT <- loadObject(pathname);
  attachLocally(dataT);
  data <- dataT;
  rm(fitSeg, dataT);
  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCN segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

force <- FALSE;

if (!exists("fitSeg") || !inherits(fitSeg, "PairedPSCBS")) {
  pathnameT <- gsub(".xdr", ",PairedPSCBS.xdr", pathname, fixed=TRUE);
  if (!force && isFile(pathnameT)) {
    verbose && enter(verbose, "Loading PairedPSCBS results");
    fitSeg <- loadObject(pathnameT);
    if (inherits(fitSeg, "PSCBS") && !inherits(fitSeg, "AbstractCBS")) {
      class(fitSeg) <- c(class(fitSeg), "AbstractCBS");
      saveObject(fitSeg, file=pathnameT);
    }
    verbose && print(verbose, fitSeg);
    verbose && exit(verbose);
  } else {
    dataT <- data.frame(CT=CT, betaT=betaT, betaN=betaN, chromosome=chromosome, x=position);
    fitSeg <- segmentByPairedPSCBS(dataT, verbose=verbose);

    verbose && enter(verbose, "Saving PairedPSCBS object");
    saveObject(fitSeg, file=pathnameT);
    verbose && exit(verbose);
  }
  rm(fitSeg2, pathnameT);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Post-segmentation of TCNs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("fitSeg2") || !inherits(fitSeg2, "PairedPSCBS")) {
  pathnameT <- gsub(".xdr", ",PairedPSCBS,psTCN.xdr", pathname, fixed=TRUE);
  if (!force && isFile(pathnameT)) {
    verbose && enter(verbose, "Loading TCN-post-segmented PairedPSCBS results");
    fitSeg2 <- loadObject(pathnameT);
    if (inherits(fitSeg2, "PSCBS") && !inherits(fitSeg2, "AbstractCBS")) {
      class(fitSeg2) <- c(class(fitSeg2), "AbstractCBS");
      saveObject(fitSeg2, file=pathnameT);
    }
    verbose && print(verbose, fitSeg2);
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Post-segmenting TCNs");
    fitSeg2 <- postsegmentTCN(fitSeg, verbose=verbose);
    verbose && print(verbose, fitSeg2);

    verbose && enter(verbose, "Saving PairedPSCBS object");
    saveObject(fitSeg2, file=pathnameT);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
  rm(pathnameT);
}

 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Calibrate in (C1,C2) space");
fitSegO <- calibrateC1C2(fitSeg2, verbose=-10);
verbose && exit(verbose);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot genomewide (C1,C2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Clim <- c(0,4);
devSet("C1C2,raw");
plotC1C2Grid(fitSeg2, Clim=Clim);
title(main="(C1,C2) - raw");

Clim <- c(-1,10);
devSet("C1C2,calibrated");
plotC1C2Grid(fitSegO, Clim=Clim);
abline(h=0); abline(v=0);
title(main="(C1,C2) - calibrated");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (C1,C2) chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (FALSE) {

devSet("C1C2,byChromosome");
par(ask=TRUE);
for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);

  fitListT <- list(
    extractByChromosome(fitSeg2, chr),
    extractByChromosome(fitSegO, chr)
  );

  Clim <- c(0,4);
  subplots(4, ncol=2, byrow=TRUE);
  par(mar=c(1,3.5,1,0.5)+1);

  for (fit in fitListT) {
    plotC1C2(fit, Clim=Clim);
    stext(side=1, pos=1, line=-1, chrTag);
    linesC1C2(fit);
  }
} # for (chr ...)

} # if (FALSE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot tracks chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
imgName <-"C1C2,tracks,byChromosome";
#devSet(imgName); par(ask=TRUE);
for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);
  imgFile <- sprintf("%s,%s.png", imgName, chrTag);
  devNew("png", imgFile, width=960, height=600);

  fitListT <- list(
    extractByChromosome(fitSeg2, chr),
    extractByChromosome(fitSegO, chr)
  );

  subplots(2, ncol=1, byrow=FALSE);
  par(mar=c(1,3.5,1,0.5)+1);
  for (fit in fitListT) {
    plotTracks(fit, tracks="tcn,c1,c2", scatter=FALSE, grid=TRUE, Clim=Clim, add=TRUE);
##    legend("topright", lwd=3, col=c("purple", "red", "blue"), c("C - Total CN", "C2 - Major CN", "C1 - Minor CN"), bty="n", horiz=TRUE);
    legend("topleft", lwd=3, col=c("purple", "red", "blue"), c("C - Total CN", "C2 - Major CN", "C1 - Minor CN"), bg="white", horiz=FALSE);
  }

  devDone();
} # for (chr ...)

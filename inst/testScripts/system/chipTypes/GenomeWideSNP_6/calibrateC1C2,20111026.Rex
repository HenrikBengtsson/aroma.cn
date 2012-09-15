library("aroma.cn");
library("PSCBS");
#sourceTo("../aroma.cn/R/PairedPSCBS.NORM.R", modifiedOnly=TRUE);


verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load ASCN data (a 10% subset)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSetName <- "GSE20584";
dataSetTags <- "ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";
dataSet <- paste(c(dataSetName, dataSetTags, "sample=0.10"), collapse=",");
chipType <- "GenomeWideSNP_6";
sampleName <- "GSM517071vGSM517072";
filename <- sprintf("%s,%s,%s.xdr", dataSetName, sampleName, dataSetTags);

rootPath <- "wholeGenomeData/";
rootPath <- Arguments$getReadablePath(rootPath);
path <- file.path(rootPath, dataSet, chipType);
path <- Arguments$getReadablePath(path);

if (!exists("data", mode="list")) {
  verbose && enter(verbose, "Loading ASCN data");
  pathname <- Arguments$getReadablePathname(filename, path=path);
  dataT <- loadObject(pathname);
  attachLocally(dataT);
  data <- dataT;
  rm(fitList, dataT);
  verbose && exit(verbose);
}


if (!exists("fitList")) {
  fitList <- list();
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCN segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

force <- FALSE;

rootPath <- "pscbsData/";
rootPath <- Arguments$getWritablePath(rootPath);

if (!is.element("raw", names(fitList))) {
  path <- file.path(rootPath, dataSet, chipType);
  path <- Arguments$getReadablePath(path, mustExists=FALSE);
  filenameT <- gsub(".xdr", ",PairedPSCBS.xdr", filename, fixed=TRUE);
  pathnameT <- file.path(path, filenameT);
  if (!force && isFile(pathnameT)) {
    verbose && enter(verbose, "Loading PairedPSCBS results");
    fit <- loadObject(pathnameT);
    verbose && print(verbose, fit);
    verbose && exit(verbose);
  } else {
    dataT <- data.frame(chromosome=chromosome, x=position, 
                        CT=CT, betaT=betaT, betaN=betaN);
    fit <- segmentByPairedPSCBS(dataT, verbose=verbose);

    verbose && enter(verbose, "Saving PairedPSCBS object");
    path <- Arguments$getWritablePath(path);
    saveObject(fit, file=pathnameT);
    verbose && exit(verbose);
  }
  rm(filenameT, pathnameT);

  fitList$raw <- fit;
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!is.element("calib", names(fitList))) {
  verbose && enter(verbose, "Calibrate in (C1,C2) space");
  fitO <- calibrateC1C2(fitList$raw, verbose=-10);
  fitList$calib <- fitO;
  verbose && exit(verbose);
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot genomewide (C1,C2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
colT <- "#66666633";
Clim <- c(0,4);
devSet("C1C2,raw");
plotC1C2Grid(fitSeg, Clim=Clim, col=colT);
#linesC1C2(fitSeg, col=colT);
title(main="(C1,C2) - raw");

Clim <- c(-1,10);
devSet("C1C2,calibrated");
plotC1C2Grid(fitSegO, Clim=Clim, col=colT);
abline(h=0); abline(v=0);
title(main="(C1,C2) - calibrated");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (C1,C2) chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (FALSE) {

devSet("C1C2,byChromosome");
#par(ask=TRUE);
for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);

  fitListT <- list(
    extractByChromosome(fitSeg, chr),
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
Clim <- c(0,4.5);
for (chr in 3) {
  chrTag <- sprintf("chr%02d", chr);
  for (key in names(fitList)) {
    fit <- extractChromosome(fitList[[key]], chr);

    if (key == "calib" && chr == 3) {
      segs <- fit$output;
      segs["18","c1Mean"] <- 0.04; segs["18","c2Mean"] <- 2.06; 
      segs["18","tcnMean"] <- segs["18","c1Mean"] + segs["18","c2Mean"];
      segs["19","c1Mean"] <- 2.015; segs["19","c2Mean"] <- 2.055;
      fit$output <- segs;
    }

    tags <- c(key, chrTag);
    toPNG(sampleName, tags=tags, width=960, aspectRatio=0.35, par=list(mar=c(1.2,0.7,1.0,0.5)+0.5, mgp=c(1.3,0.5,0), las=1, cex=2), {
      plotTracks(fit, tracks="tcn,c1,c2", scatter=c(), calls=c(), lwd=6, grid=TRUE, Clim=Clim, ylab="");
      drawLevels(fit, what="tcn", col="purple", lwd=6);
      drawLevels(fit, what="c1", col="blue", lwd=6);
      drawLevels(fit, what="c2", col="red", lwd=6);
#     legend("topleft", lwd=6, col=c("purple", "red", "blue"), c("C - Total CN", "C2 - Major CN", "C1 - Minor CN"), bg="white", horiz=FALSE, bty="n");
      stext(side=3, pos=0, sprintf("%s:%s", dataSetName, sampleName));
    });
  } # for (key ...)
} # for (chr ...)

graphics.off();

###########################################################################
# HISTORY:
# 2011-10-26
# o Updated to work with the released version of PSCBS. Previously
#   it was based on a pre-release version of PSCBS.
# [...]
# 2010-11-11
# o Cool calibration plots for GSE20584,GSM517071vGSM517072 on Chr3.
# 2010-xx-xx
# o Created.
###########################################################################

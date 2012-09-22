library("aroma.cn");
library("PSCBS");
library("R.devices");
library("R.menu");
verbose <- Arguments$getVerbose(-10);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("doPlots", "PairedPSCBS", function(fit, sampleName=NULL, tags=NULL, ...) {
  # Argument 'sampleName':
  if (is.null(sampleName)) {
    sampleName <- sampleName(fit);
  }
  stopifnot(!is.null(sampleName));

  nCPsTag <- sprintf("#CPs=%d", nbrOfChangePoints(fit));
  toPNG(sampleName, tags=c("(C1,C2)", nCPsTag, tags), width=800, {
    plotC1C2Grid(fit);
    linesC1C2(fit);
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, nCPsTag);
    stext(side=4, pos=0, dataSet, cex=0.7);
    stext(side=4, pos=1, chipType, cex=0.7);
  });
  
  
  toPNG(sampleName, tags=c("tracks", nCPsTag, tags), width=1200, aspectRatio=0.25, {
    plotTracks(fit, tracks="tcn,c1,c2");
    stext(side=4, pos=0, sampleName);
    stext(side=4, pos=1, nCPsTag);
  });
}) # doPlots()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup Paired PSCBS segmentation data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "pscbsData";
path <- Arguments$getReadablePath(rootPath);

dataSets <- list.files(rootPath);
if (length(dataSets) > 1) {
 dataSet <- textMenu(dataSets, value=TRUE);
} else {
 dataSet <- dataSets[1];
}

path <- file.path(rootPath, dataSet);
path <- Arguments$getReadablePath(path);
chipTypes <- list.files(path);
if (length(chipTypes) > 1) {
 chipType <- textMenu(chipTypes, value=TRUE);
} else {
 chipType <- chipTypes[1];
}

ds <- PairedPSCBSFileSet$byName(dataSet, chipType=chipType);
print(ds);
dsName <- getName(ds);

if (length(ds) == 0) {
 throw("No PairedPSCBS data file found.")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (length(ds) > 1) {
 ii <- textMenu(getNames(ds));
} else {
 ii <- 1L;
}

df <- getFile(ds, ii);
fit <- loadObject(df);
sampleName <- getName(df);

fit0 <- fit;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Configure report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
figPath <- file.path("figures", dataSet);
options("devEval/args/path"=figPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (C1,C2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
doPlots(fit);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prune change points using dynamic programming
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("segList", mode="list")) {
  segList <- seqOfSegmentsByDP(fit, verbose=-10);
  modelFit <- attr(segList, "modelFit");
  modelFit$seqOfSegmentsByDP <- NULL;
  str(modelFit);
}


toPNG(dataSet, tags=c(sampleName, "DP", "RSEvsCPs"), width=800, aspectRatio=0.7, {
  plot(modelFit$nbrOfChangePoints, modelFit$rse,
       xlab="Number of change points", ylab="RSE");
  stext(side=3, pos=0, sampleName);
  stext(side=4, pos=0, dataSet, cex=0.7);
  stext(side=4, pos=1, chipType, cex=0.7);
});


nbrOfCPs <- c(100, 50, 25)[1:2];
if (!exists("fitList", mode="list")) {
  fitList <- list();
}
for (kk in seq(along=nbrOfCPs)) {
  key <- sprintf("nbrOfCPs=%d", nbrOfCPs[kk]);

  verbose && enter(verbose, sprintf("Change point set #%d ('%s') of %d", kk, key, length(nbrOfCPs)));

  verbose && cat(verbose, "Number of change points: ", nbrOfCPs[kk]);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Pruning CPs via dynamic programming
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- fitList[[key]];
  if (is.null(fitT)) {
    verbose && enter(verbose, "Resegmenting");
    knownSegments <- segList[[nbrOfCPs[kk]+1L]];
    fitT <- resegment(fit, knownSegments=knownSegments, undoTCN=+Inf, undoDH=+Inf);
    fitList[[key]] <- fitT;
    verbose && exit(verbose);
  }
  sampleName(fitT) <- sampleName(fit);
  fitDP <- fitT;
  doPlots(fitDP);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Deshear
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fast naive calling of AB
##   fitAB <- fitDP;
##   deltaAB <- estimateDeltaAB(fitAB);
##   segs <- fitAB$output;
##   abCall <- segs$dhMean <= deltaAB;
##   segs$dhMean[abCall] <- 0;
##   segs$abCall <- abCall;
##   fitAB$output <- segs;
  fitD <- deShearC1C2(fitDP);
  doPlots(fitD, tags="deShear");

  verbose && exit(verbose);
} # for (kk ...)

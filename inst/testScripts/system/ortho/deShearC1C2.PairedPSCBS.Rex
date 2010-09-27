library("aroma.cn");
library("psCBS");
#sourceTo("../aroma.cn/R/PairedPSCBS.NORM.R", modifiedOnly=TRUE);

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load ASCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "wholeGenomeData";
dataSet <- "GSE20584,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";
chipType <- "GenomeWideSNP_6";
path <- file.path(rootPath, dataSet, chipType);
filename <- "GSE20584,GSM517071vGSM517072,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY.xdr";
pathname <- Arguments$getReadablePathname(filename, path=path);

if (!exists("data", mode="list")) {
  verbose && enter(verbose, "Loading ASCN data");
  data <- loadObject(pathname);
  rm(fitList);
  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCN segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathnameT <- gsub(".xdr", ",fitList.xdr", pathname, fixed=TRUE);
if (!exists("fitList", mode="list")) {
  if (isFile(pathnameT)) {
    verbose && enter(verbose, "Loading PSCN results");
    fitList <- loadObject(pathnameT);
    verbose && exit(verbose);
  } else {
    fitList <- list();
  }
}

wasUpdated <- FALSE;
for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);
  verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", chr, chrTag, 22));

  fit <- fitList[[chrTag]];
  if (is.null(fit)) {
    verbose && enter(verbose, "Segmenting");
    units <- whichVector(data$chromosome == chr);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);
    dataC <- data[units,,drop=FALSE];
    attachLocally(dataC);
    x <- position;
    # AD HOC: Robustification
    CT[CT < 0] <- 0;
    CT[CT > 30] <- 30;
    fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, 
                                chromosome=chr, x=x, verbose=-10);
    fitList[[chrTag]] <- fit;
    rm(units, CT, betaT, betaN, chr, x, dataC);
    wasUpdated <- TRUE;
    verbose && exit(verbose);
  } else {
    verbose && cat(verbose, "Already segmented.");
  }
  rm(fit);  
  verbose && exit(verbose);
} # for (chr ...)

if (!isFile(pathnameT) || wasUpdated) {
  verbose && enter(verbose, "Saving PSCN results");
  saveObject(fitList, file=pathnameT);
  rm(fitSeg);
  verbose && exit(verbose);
}
rm(pathnameT);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reduce list
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("fitSeg") || !inherits(fitSeg, "PairedPSCBS")) {
  pathnameT <- gsub(".xdr", ",PairedPSCBS.xdr", pathname, fixed=TRUE);
  if (isFile(pathnameT)) {
    verbose && enter(verbose, "Loading PairedPSCBS results");
    fitSeg <- loadObject(pathnameT);
    verbose && print(verbose, fitSeg);
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Reducing list of PSCN results to PairedPSCBS object");
    fitSeg <- Reduce(append, fitList);
    verbose && print(verbose, fitSeg);

    verbose && enter(verbose, "Saving PairedPSCBS object");
    saveObject(fitSeg, file=pathnameT);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
  rm(fitSeg2);
  rm(pathnameT);
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Post-segmentation of TCNs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("fitSeg2") || !inherits(fitSeg2, "PairedPSCBS")) {
  pathnameT <- gsub(".xdr", ",PairedPSCBS,psTCN.xdr", pathname, fixed=TRUE);
  if (isFile(pathnameT)) {
    verbose && enter(verbose, "Loading TCN-post-segmented PairedPSCBS results");
    fitSeg2 <- loadObject(pathnameT);
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
# Normalize BAF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("fitSegN") || !inherits(fitSegN, "PairedPSCBS")) {
  pathnameT <- gsub(".xdr", ",PairedPSCBS,psTCN,nBAF.xdr", pathname, fixed=TRUE);
  if (isFile(pathnameT)) {
    verbose && enter(verbose, "Loading BAF-normalized TCN-post-segmented PairedPSCBS results");
    fitSegN <- loadObject(pathnameT);
    verbose && print(verbose, fitSegN);
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Normalizing BAF by regions");
    fitSegN <- normalizeBAFsByRegions(fitSeg2, verbose=-10);
    verbose && print(verbose, fitSegN);

    verbose && enter(verbose, "Saving PairedPSCBS object");
    saveObject(fitSegN, file=pathnameT);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
  rm(pathnameT);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call allelic balance (AB) and adjust for AB biases
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("fitSegC") || !inherits(fitSegC, "PairedPSCBS")) {
  pathnameT <- gsub(".xdr", ",PairedPSCBS,psTCN,nBAF,AB.xdr", pathname, fixed=TRUE);
  if (isFile(pathnameT)) {
    verbose && enter(verbose, "Loading AB-called BAF-normalized TCN-post-segmented PairedPSCBS results");
    fitSegC <- loadObject(pathnameT);
    verbose && print(verbose, fitSegC);
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Normalizing for allelic-balance biases");
    fitSegC <- callAllelicBalanceByBAFs(fitSegN, verbose=-10);
    ww <- which(fitSegC$output$ab.call);
    fitSegC$output[ww, "dh.mean"] <- 0;
    verbose && print(verbose, fitSegC);

    verbose && enter(verbose, "Saving PairedPSCBS object");
    saveObject(fitSegC, file=pathnameT);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
}
 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Deshear
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Deshearing in (C1,C2) space");
fitSegS <- deShearC1C2(fitSegC, verbose=-10);
verbose && exit(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot genomewide (C1,C2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2");
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);
plotC1C2(fitSeg, Clim=Clim);
title(main="(C1,C2) - raw");
plotC1C2(fitSegS, Clim=Clim);
title(main="(C1,C2) - desheared");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (C1,C2) chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);

  fit <- extractByChromosome(fitSeg, chr);
  fit2 <- extractByChromosome(fitSeg2, chr);
  fitS <- extractByChromosome(fitSegS, chr);

  devSet(sprintf("C1C2,%s", chrTag));
  Clim <- c(0,4);
  subplots(4, ncol=2, byrow=TRUE);
  par(mar=c(1,3.5,1,0.5)+1);
  plotC1C2(fit, Clim=Clim);
  stext(side=1, pos=1, line=-1, chrTag);
  linesC1C2(fit);
  title(main="(C1,C2)");

  plotC1C2(fitS, Clim=Clim);
  stext(side=1, pos=1, line=-1, chrTag);
  linesC1C2(fitS);
  title(main="(C1,C2) - desheared");
} # for (chr ...)


devSet("tracks");
subplots(10, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plot(fit, add=TRUE)
plot(fitS, add=TRUE)



stop()
stop()
stop()



if (!exists("fitNList", mode="list")) {
  fitNList <- list();
}

for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);
  fitN <- fitNList[[chrTag]];
  if (is.null(fitN)) {
    fit2 <- fit2List[[chrTag]];
    fitN <- normalizeBAFsByRegions(fit2, verbose=-10);
    fitNList[[chrTag]] <- fitN;
  }
  rm(fitN);  
} # for (chr ...)


c1c2 <- extractC1C2(fitList);
c1c2N <- extractC1C2(fitNList);

fit2D <- fitLoess2D(X=c1c2, Y=c1c2N);
c1c2N2 <- normalizeLoess2D(c1c2, fit=fit2D);

AB <- C1C2toAB(c1c2N2);

devSet("C1C2");
plot(c1c2N2);
lines(c1c2N2);
plot(AB[,c("l","b")], cex=sqrt(c1c2N2[,"w"])+1/2);


fitT <- fitList[[1]];
C1C2 <- extractC1C2(fitT);

# Orthogonalize
devSet("ortho");
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);

stop()

devSet("tracks");
subplots(10, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plot(fitT, add=TRUE);
plot(fitO, add=TRUE);

plotC1C2(fitT, col="gray", Clim=Clim);
linesC1C2(fitT, col="gray");
pointsC1C2(fitO);
linesC1C2(fitO);
title(main="(C1,C2) - both");

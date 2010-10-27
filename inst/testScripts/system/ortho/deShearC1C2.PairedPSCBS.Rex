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
  rm(fitSeg);
  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCN segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

if (!exists("fitSeg") || !inherits(fitSeg, "PairedPSCBS")) {
  pathnameT <- gsub(".xdr", ",PairedPSCBS.xdr", pathname, fixed=TRUE);
  if (isFile(pathnameT)) {
    verbose && enter(verbose, "Loading PairedPSCBS results");
    fitSeg <- loadObject(pathnameT);
    verbose && print(verbose, fitSeg);
    verbose && exit(verbose);
  } else {
    fitSeg <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, chromosome=chromosome, x=position, verbose=verbose);

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

fitSegSS <- deShearC1C2(fitSegS, dirs="X");
fitSegSS <- deShearC1C2(fitSegSS, dirs="|,-,X");

dC2 <- estimateC2Bias(fitSegSS);
fitSegSSS <- translateC1C2(fitSegSS, dC2=-dC2);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot genomewide (C1,C2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2");
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);
plotC1C2(fitSeg2, Clim=Clim);
title(main="(C1,C2)");
plotC1C2(fitSegN, Clim=Clim);
title(main="(C1,C2) - BAF normalized");
plotC1C2(fitSegS, Clim=Clim);
title(main="(C1,C2) - desheared");
plotC1C2(fitSegSS, Clim=Clim);
title(main="(C1,C2) - desheared + rescaled");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot genomewide (C1,C2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2,grid");
plotC1C2Grid(fitSegS, main="(C1,C2) - desheared");

devSet("C1C2,grid,scaled");
plotC1C2Grid(fitSegSS, main="(C1,C2) - desheared + rescaled");

devSet("C1C2,grid,scaled");
plotC1C2Grid(fitSegSSS, main="(C1,C2) - desheared + rescaled + shift");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (C1,C2) chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2,byChromosome");
par(ask=TRUE);
for (chr in 1) {
  chrTag <- sprintf("chr%02d", chr);

  fitListT <- list(
    extractByChromosome(fitSeg2, chr),
    extractByChromosome(fitSegN, chr),
    extractByChromosome(fitSegS, chr),
    extractByChromosome(fitSegSS, chr),
    extractByChromosome(fitSegSSS, chr)
  );

  Clim <- c(0,4);
  subplots(9, ncol=3, byrow=TRUE);
  par(mar=c(1,3.5,1,0.5)+1);

  for (fit in fitListT) {
    plotC1C2(fit, Clim=Clim);
    stext(side=1, pos=1, line=-1, chrTag);
    linesC1C2(fit);
  }
} # for (chr ...)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (x,C1,C2) chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2,byChromosome,tracks");
par(ask=TRUE);
for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);

  fitListT <- list(
    extractByChromosome(fitSeg2, chr),
    extractByChromosome(fitSegN, chr),
    extractByChromosome(fitSegS, chr),
    extractByChromosome(fitSegSS, chr),
    extractByChromosome(fitSegSSS, chr)
  );

  c1c2List <- lapply(fitListT, FUN=function(fit) {
    xx <- as.data.frame(fit)[,c("tcn.loc.start", "tcn.loc.end"),drop=FALSE];
    c1c2 <- extractC1C2(fit)[,1:2];
    cbind(xx, c1c2);
  });

  c1c2 <- c1c2List[[1]];
  xrange <- as.matrix(c1c2[,c("tcn.loc.start", "tcn.loc.end")]);
  xrange <- range(xrange, na.rm=TRUE);

  subplots(length(c1c2List), ncol=1, byrow=FALSE);
  par(mar=c(1,3.5,1,0.5)+1);

  for (kk in seq(along=c1c2List)) {
    c1c2 <- c1c2List[[kk]];
    plot(NA, xlim=xrange, ylim=Clim);
    stext(side=3, pos=1, chrTag);
    dd <- apply(c1c2, MARGIN=1, FUN=function(df) {
      xx <- as.numeric(df[1:2]);
      c1 <- rep(as.numeric(df[3]), times=2);
      c2 <- rep(as.numeric(df[4]), times=2);
      c <- c1+c2;
      lines(xx, c, lwd=2, col=1);
      lines(xx, c1, lwd=2, col=2);
      lines(xx, c2, lwd=2, col=3);
    });
  }
} # for (chr ...)

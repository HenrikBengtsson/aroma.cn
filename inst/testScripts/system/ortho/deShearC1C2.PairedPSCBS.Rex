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
plotC1C2(fitSeg2, Clim=Clim);
title(main="(C1,C2)");
plotC1C2(fitSegN, Clim=Clim);
title(main="(C1,C2) - BAF normalized");
plotC1C2(fitSegC, Clim=Clim);
title(main="(C1,C2) - AB adjusted");
plotC1C2(fitSegS, Clim=Clim);
title(main="(C1,C2) - desheared");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot genomewide (C1,C2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2,grid");
Clim <- c(0,4);
fit <- fitSegS;
plotC1C2(fit, Clim=Clim);
title(main="(C1,C2) - desheared");
data <- extractC1C2(fit);
n <- data[,4, drop=TRUE];
n <- sqrt(n);
w <- n/sum(n, na.rm=TRUE);
adjust <- 0.2;
for (cc in 1:2) {
  y <- data[,cc];
  ok <- is.finite(y) & is.finite(w);
  y <- y[ok];
  wt <- w[ok]/sum(w[ok]);
  d <- density(y, weights=wt, adjust=adjust);
  draw(d, side=cc, height=0.3, col="gray", lwd=2, xpd=FALSE);
  p <- findPeaksAndValleys(y, weights=wt, adjust=adjust);
  p <- subset(p, type == "peak" & density > 0.05);
  p <- p[order(p$density, decreasing=TRUE),,drop=FALSE];
  p <- head(p, n=8);
  if (cc == 1) {
    abline(v=p$x, lty=3, col="gray");
  } else {
    abline(h=p$x, lty=3, col="gray");
  }
}
box();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (C1,C2) chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2,byChromosome");
par(ask=TRUE);
for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);

  fitList <- list(
    extractByChromosome(fitSeg2, chr),
    extractByChromosome(fitSegN, chr),
    extractByChromosome(fitSegC, chr),
    extractByChromosome(fitSegS, chr)
  );

  Clim <- c(0,4);
  subplots(4, ncol=2, byrow=TRUE);
  par(mar=c(1,3.5,1,0.5)+1);

  for (fit in fitList) {
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

  fitList <- list(
    extractByChromosome(fitSeg2, chr),
    extractByChromosome(fitSegN, chr),
    extractByChromosome(fitSegC, chr),
    extractByChromosome(fitSegS, chr)
  );

  c1c2List <- lapply(fitList, FUN=function(fit) {
    xx <- as.data.frame(fit)[,c("tcn.loc.start", "tcn.loc.end"),drop=FALSE];
    c1c2 <- extractC1C2(fit)[,1:2];
    cbind(xx, c1c2);
  });

  c1c2 <- c1c2List[[1]];
  xrange <- as.matrix(c1c2[,c("tcn.loc.start", "tcn.loc.end")]);
  xrange <- range(xrange, na.rm=TRUE);

  subplots(length(c1c2List), ncol=1, byrow=FALSE);
  par(mar=c(1,3.5,1,0.5)+1);

  for (c1c2 in c1c2List) {
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

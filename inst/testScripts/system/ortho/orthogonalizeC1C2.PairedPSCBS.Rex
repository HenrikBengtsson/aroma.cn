library("aroma.cn");
library("psCBS");
#sourceTo("../aroma.cn/R/PairedPSCBS.NORM.R", modifiedOnly=TRUE);

pathname <- "GSE20584,GSM517071vGSM517072,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY.xdr";
pathname2 <- gsub(".xdr", ",fitList.xdr", pathname, fixed=TRUE);

if (!exists("data", mode="list")) {
  # Load example ASCN data
  data <- loadObject(pathname);
  rm(dataList);
}

if (!exists("dataList", mode="list")) {
  dataList <- list();
  rm(fitList);
}

if (!exists("fitList", mode="list")) {
  if (isFile(pathname2)) {
    fitList <- loadObject(pathname2);
  } else {
    fitList <- list();
  }
}

if (!exists("fit2List", mode="list")) {
  fit2List <- list();
}

if (!exists("fitNList", mode="list")) {
  fitNList <- list();
}

for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);

  dataC <- dataList[[chrTag]];
  if (is.null(dataC)) {
    units <- whichVector(data$chromosome == chr);
    str(units);
    dataC <- data[units,,drop=FALSE];
    dataList[[chrTag]] <- dataC;
    rm(units);
  }

  fit <- fitList[[chrTag]];
  if (is.null(fit)) {
    attachLocally(dataC);
    x <- position;
    # AD HOC: Robustification
    CT[CT < 0] <- 0;
    CT[CT > 30] <- 30;
    fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, 
                                chromosome=chr, x=x, verbose=-10);
    fitList[[chrTag]] <- fit;
  }
  rm(dataC, fit);  
} # for (chr ...)


for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);
  fit2 <- fit2List[[chrTag]];
  if (is.null(fit2)) {
    fit <- fitList[[chrTag]];
    fit2 <- postsegmentTCN(fit);
    fit2List[[chrTag]] <- fit2;
  }
  rm(fit2);  
} # for (chr ...)


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


if (!isFile(pathname2)) {
  saveObject(fitList, file=pathname2);
}

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
fitO <- orthogonalizeC1C2(fitT, debugPlot=TRUE, verbose=-10);

stop()

devSet("tracks");
subplots(10, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plot(fitT, add=TRUE);
plot(fitO, add=TRUE);

devSet("C1C2");
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);
plotC1C2(fitT, Clim=Clim);
linesC1C2(fitT);
title(main="(C1,C2) - adjusted");

plotC1C2(fitO, Clim=Clim);
linesC1C2(fitO);
title(main="(C1,C2) - orthogonalized");

plotC1C2(fitT, col="gray", Clim=Clim);
linesC1C2(fitT, col="gray");
pointsC1C2(fitO);
linesC1C2(fitO);
title(main="(C1,C2) - both");


library("aroma.cn");
library("PSCBS");
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

if (!exists("fitNList", mode="list")) {
  fitNList <- list();
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation for each chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  fitN <- fitNList[[chrTag]];
  if (is.null(fitN)) {
    fitN <- normalizeBAFsByRegions(fit, verbose=-10);
    fitNList[[chrTag]] <- fitN;
  }

  rm(dataC, fit, fitN);  
} # for (chr ...)

# Save results?
if (!isFile(pathname2)) {
  saveObject(fitList, file=pathname2);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit (C1,C2) normalization function from DH vs TCN effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c1c2 <- extractC1C2(fitList);
c1c2N <- extractC1C2(fitNList);
fit2D <- fitLoess2D(X=c1c2, Y=c1c2N);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize (C1,C2) accordingly
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c1c2N2 <- normalizeLoess2D(c1c2, fit=fit2D);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot normalization results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2");
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(3.5,3.5,0.5,0.5)+1);
plot(c1c2); title(main="(C1,C2)");
plot(c1c2N); title(main="(C1,C2) - adjusted");
plot(c1c2N2); title(main="(C1,C2) - normalized");

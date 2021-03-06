if (interactive()) savehistory();
library("aroma.cn");
library(R.menu);

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

## dataSet <- "TCGA,GBM,onePair";
dataSet <- "TCGA,OV,testSet,pairs,Broad,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";

if (interactive()) {
  dataSet <- textMenu(dataSets, value=TRUE);
} else {
  ## dataSet <- "TCGA,GBM,onePair";
  dataSet <- "TCGA,OV,testSet,pairs,Broad,ismpolish";
  dataSet <- "TCGA,OV,testSet,pairs,Stanford";
  dataSet <- "TCGA,OV,testSet,pairs,Broad,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";
}
print(dataSet);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the raw (tumor,normal) data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ds <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
setFullNamesTranslator(ds, function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
});
print(ds);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract the normals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
types <- sapply(ds, FUN=function(df) getTags(df)[1]);
normals <- grep("(10|11)[A-Z]", types);
dsN <- extract(ds, normals);
print(dsN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Naive genotype calling algorithm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
rootPath <- "totalAndFracBData";
fullname <- paste(c(getFullName(dsN), "NGC"), collapse=",");
chipType <- getChipType(dsN, fullname=FALSE);
outPath <- file.path("callData", fullname, chipType);

units <- NULL;
if (is.null(units)) {
  df <- getFile(dsN, 1);
  units <- seq_len(nbrOfUnits(df));
  rm(df);
}

adjust <- 1.5;
type <- NULL; rm(type); # To please R CMD check

# Identify units on ChrX and ChrY
ugp <- getAromaUgpFile(dsN);
units23 <- getUnitsOnChromosome(ugp, 23);
is23 <- is.element(units, units23);
units24 <- getUnitsOnChromosome(ugp, 24);
is24 <- is.element(units, units24);

for (kk in seq(dsN)) {
  dfN <- getFile(dsN, kk);

  tags <- getTags(dfN);
  tags <- setdiff(tags, "fracB");
  tags <- c(tags, "genotypes");
  fullname <- paste(c(getName(dfN), tags), collapse=",");
  filename <- sprintf("%s.acf", fullname);
  gcPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE);

  csTags <- c(tags, "confidenceScores");
  fullname <- paste(c(getName(dfN), csTags), collapse=",");
  filename <- sprintf("%s.acf", fullname);
  csPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE);
  
  if (isFile(gcPathname) && isFile(csPathname)) {
    next;
  }

  betaN <- dfN[units,1,drop=TRUE];

  # Call gender
  gender <- callXXorXY(betaN[is23], betaN[is24], adjust=adjust, from=0, to=1);

  # Call genotypes 
  naValue <- as.double(NA);
  fit <- NULL;
  mu <- rep(naValue, times=length(units));
  cs <- rep(naValue, times=length(units));

  if (gender == "XY") {
    # All but ChrX & ChrY in male
    isDiploid <- (!(is23 | is24));
    use <- which(isDiploid);
    muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, from=0, to=1,
                                          verbose=less(verbose,10));
    fit <- attr(muT, 'modelFit');
    mu[use] <- muT;
    use <- which(!isDiploid);
    muT <- callNaiveGenotypes(betaN[use], cn=1, adjust=adjust, from=0, to=1,
                                           verbose=less(verbose,10));
    mu[use] <- muT;
  } else {
    # All but ChrY in female
    isDiploid <- (!is24);
    use <- which(isDiploid);
    muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, from=0, to=1,
                                          verbose=less(verbose,10));
    fit <- attr(muT, 'modelFit');
    mu[use] <- muT;
  }
  print(table(mu, exclude=NULL));

  # Translate genotype calls in fracB space to (AA,AB,BB,...)
  calls <- rep(as.character(NA), times=length(mu));
  calls[mu ==   0] <- "AA";
  calls[mu == 1/2] <- "AB";
  calls[mu ==   1] <- "BB";
  print(table(calls, exclude=NULL));

  # Calculate confidence scores
  a <- fit$x[1];
  b <- fit$x[2]; 
  cs[isDiploid] <- rowMins(abs(cbind(betaN[isDiploid]-a, betaN[isDiploid]-b)));
  print(table(mu, exclude=NULL));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Writing genotype calls (via temporary file)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pathname <- gcPathname;
  pathnameT <- sprintf("%s.tmp", pathname);
  nbrOfUnits <- nbrOfUnits(dfN);
  gfN <- AromaUnitGenotypeCallFile$allocate(pathnameT, platform=getPlatform(dfN), chipType=getChipType(dfN), nbrOfRows=nbrOfUnits);
  footer <- readFooter(gfN);
  footer$method <- "NaiveGenotypeCaller";
  writeFooter(gfN, footer);
  rm(footer);

  updateGenotypes(gfN, units=units, calls=calls);
  rm(calls);

  res <- file.rename(pathnameT, pathname);
  if (!isFile(pathname)) {
    throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
  }
  if (isFile(pathnameT)) {
    throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
  }
  rm(pathnameT);

  gfN <- AromaUnitGenotypeCallFile(pathname);
  print(gfN);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Writing confidence scores (via temporary file)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- csPathname;
  pathnameT <- sprintf("%s.tmp", pathname);
  nbrOfUnits <- nbrOfUnits(dfN);
  csfN <- AromaUnitSignalBinaryFile$allocate(pathnameT, platform=getPlatform(dfN), chipType=getChipType(dfN), nbrOfRows=nbrOfUnits, types="double", size=4, signed=TRUE);
  footer <- readFooter(csfN);
  footer$method <- "NaiveGenotypeConfidenceScoreEstimator";
  writeFooter(csfN, footer);
  rm(footer);

  csfN[units, 1] <- cs
  rm(cs);

  res <- file.rename(pathnameT, pathname);
  if (!isFile(pathname)) {
    throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
  }
  if (isFile(pathnameT)) {
    throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
  }
  rm(pathnameT);

  rm(units);

} # for (kk ...)

gcN <- AromaUnitGenotypeCallSet$byName(dataSet, tags="NGC", chipType="*");
print(gcN);

csN <- AromaUnitSignalBinarySet$byName(dataSet, tags="NGC", chipType="*", pattern="confidenceScores", paths="callData");
print(csN);
############################################################################
## HISTORY:
## 2010-07-25
## o Now use callNaiveGenotypes (wrapper) to call naive genotypes.
## 2009-11-05
## o BUG FIX: Now use findPeakAndValleys(..., from=0, to=1) to prevent
## points outside [0,1] to fool density estimation.
## 2009-10-28
## o Now uses R.menu to choose data set.
## 2009-06-20
## o Added confidence scores.
## o Created from test20090429a,NGC.Rex.
############################################################################

if (interactive()) savehistory();
library("aroma.cn");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";

## dataSet <- "TCGA,GBM,onePair";
dataSet <- "TCGA,OV,testSet,pairs,Broad";

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
  units <- seq(length=nbrOfUnits(df));
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
  pathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE);
  if (isFile(pathname)) {
    next;
  }

  betaN <- dfN[units,1,drop=TRUE];
  keep <- is.finite(betaN);
  betaN <- betaN[keep];
  unitsN <- units[keep];
  is23k <- is23[keep];
  is24k <- is24[keep];

  # Infer male or female
  fit <- findPeaksAndValleys(betaN[is23k], adjust=adjust);
  isXY <- (sum(fit$type == "peak") == 2);
  if (sum(is24k) > 100) {
    fit <- findPeaksAndValleys(betaN[is24k], adjust=adjust);
    is24XY <- (sum(fit$type == "peak") == 2);
    if (is24XY != isXY) {
      throw("ChrX and ChrY is inconsistent");
    }
  }

  naValue <- as.double(NA);
  mu <- rep(naValue, times=length(unitsN));

  # Call genotypes for autosomal+ChrX chromosomes
  if (isXY) {
    # All but ChrX & ChrY in male
    isDiploid <- (!(is23k | is24k));
  } else {
    # All but ChrY in female
    isDiploid <- (!is24k);
  }
  fit <- findPeaksAndValleys(betaN[isDiploid], adjust=adjust);
  fit <- subset(fit, type == "valley");
  # Sanity check
  nbrOfGenotypeGroups <- nrow(fit)+1; 
  stopifnot(nbrOfGenotypeGroups == 3);
  a <- fit$x[1];
  b <- fit$x[2]; 
  mu[isDiploid & (betaN < a)] <- 0;
  mu[isDiploid & (betaN > b)] <- 1; 
  mu[isDiploid & (a <= betaN & betaN <= b)] <- 1/2;
  print(table(mu, exclude=NULL));

  # Call genotypes for ChrX+ChrY chromosomes
  if (isXY) {
    fit <- findPeaksAndValleys(betaN[!isDiploid], adjust=adjust);
    fit <- subset(fit, type == "valley");
    nbrOfGenotypeGroups <- nrow(fit)+1; 
    stopifnot(nbrOfGenotypeGroups == 2);
    a <- fit$x[1];
    mu[!isDiploid & (betaN <= a)] <- 0;
    mu[!isDiploid & (betaN  > a)] <- 1;
  } else {
    # ChrY on female
    mu[!isDiploid] <- NA;
  }
  print(table(mu, exclude=NULL));

  # Translate genotype calls in fracB space to (AA,AB,BB,...)
  calls <- rep(as.character(NA), times=length(mu));
  calls[mu ==   0] <- "AA";
  calls[mu == 1/2] <- "AB";
  calls[mu ==   1] <- "BB";
  print(table(calls, exclude=NULL));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Writing genotype calls (via temporary file)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pathnameT <- sprintf("%s.tmp", pathname);
  nbrOfUnits <- nbrOfUnits(dfN);
  gfN <- AromaUnitGenotypeCallFile$allocate(pathnameT, platform=getPlatform(dfN), chipType=getChipType(dfN), nbrOfRows=nbrOfUnits);
  footer <- readFooter(gfN);
  footer$method <- "NaiveGenotypeCaller";
  writeFooter(gfN, footer);
  rm(footer);

  updateGenotypes(gfN, units=unitsN, calls=calls);
  rm(unitsN, calls);

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
} # for (kk ...)

gcN <- AromaUnitGenotypeCallSet$byName(dataSet, tags="NGC", chipType="*");
print(gcN);

############################################################################
## HISTORY:
## 2009-06-09
## o BUG FIX: 'is23' and 'is24' were modified within the 'kk' loop, causing an
##   error when kk>1.
## 2009-06-08
## o Now uses 'AromaUnitGenotypeCallFile' and 'AromaUnitGenotypeCallSet'.
## 2009-04-29
## o Created.
############################################################################

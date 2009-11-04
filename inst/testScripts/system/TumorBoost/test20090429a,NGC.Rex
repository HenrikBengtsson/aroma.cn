if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

pattern <- "^TCGA,OV,testSet,pairs";
dataSets <- list.files(path=rootPath, pattern=pattern);
exclPattern <- ",TBN";
keep <- (regexpr(exclPattern, dataSets) == -1);
dataSets <- dataSets[keep];
# Sanity check
stopifnot(length(dataSets) > 0);

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

  # Call gender
  gender <- callXXorXY(betaN[is23], betaN[is24], adjust=adjust);

  # Call genotypes 
  naValue <- as.double(NA);
  mu <- rep(naValue, times=length(units));
  if (gender == "XY") {
    # All but ChrX & ChrY in male
    isDiploid <- (!(is23 | is24));
    use <- which(isDiploid);
    muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, 
                                          verbose=less(verbose,10));
    mu[use] <- muT;
    use <- which(!isDiploid);
    muT <- callNaiveGenotypes(betaN[use], cn=1, adjust=adjust, 
                                           verbose=less(verbose,10));
    mu[use] <- muT;
  } else {
    # All but ChrY in female
    isDiploid <- (!is24);
    use <- which(isDiploid);
    muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, 
                                          verbose=less(verbose,10));
    mu[use] <- muT;
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
} # for (kk ...)

gcN <- AromaUnitGenotypeCallSet$byName(dataSet, tags="NGC", chipType="*");
print(gcN);



############################################################################
# HISTORY:
# 2009-11-03
# o Updated to make use of new/updated callXXorXY() & callNaiveGenotypes().
# 2009-06-09
# o BUG FIX: 'is23' and 'is24' were modified within the 'kk' loop, 
#   causing an error when kk>1.
# 2009-06-08
# o Now uses 'AromaUnitGenotypeCallFile' and 'AromaUnitGenotypeCallSet'.
# 2009-04-29
# o Created.
############################################################################

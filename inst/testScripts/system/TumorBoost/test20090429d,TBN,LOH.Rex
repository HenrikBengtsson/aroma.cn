if (interactive()) savehistory();
library("aroma.cn");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";

dataSet <- "TCGA,GBM,onePair";
sampleName <- "TCGA-12-0620";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the raw (tumor,normal) data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ds <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
setFullNamesTranslator(ds, function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
});
print(ds);

pair <- indexOf(ds, sampleName);
stopifnot(length(pair) == 2);

# Order as (tumor,normal)
types <- sapply(extract(ds,pair), FUN=function(df) getTags(df)[1]);
o <- order(types);
types <- types[o];
pair <- pair[o];

# Extract (tumor, normal) pair
dsPair <- extract(ds, pair);
dsT <- extract(dsPair, 1);
dsN <- extract(dsPair, 2);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the normalized tumor data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsTN <- AromaUnitFracBCnBinarySet$byName(dataSet, tags="TBN", chipType="*", paths=rootPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the genotype call set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
gsN <- AromaUnitFracBCnBinarySet$byName(dataSet, tags="Birdseed", chipType="*", paths=rootPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create an list of matched data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsList <- list(normal=dsN, tumor=dsT, tumorN=dsTN, callsN=gsN);
dsList <- lapply(dsList, FUN=function(ds) {
  ds <- setFullNamesTranslator(ds, function(names, ...) {
    pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
    gsub(pattern, "\\1,\\2,\\3", names);
  });
  idxs <- indexOf(ds, getNames(dsList$normal));
  extract(ds, idxs);
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- getAromaUgpFile(dsList$tumor);
chromosome <- 1;
units <- getUnitsOnChromosome(ugp, chromosome=chromosome);
pos <- getPositions(ugp, units=units);


kk <- 1;
dfList <- lapply(dsList, FUN=getFile, kk);
idx <- grep("callsN", names(dfList));
cf <- dfList[[idx]];
dfList <- dfList[-idx];
callsN <- cf[units,1,drop=TRUE];

# Stratify on heterozygote normal genotypes
isHet <- whichVector(callsN == 1/2);

fracBList <- lapply(dfList, FUN=function(df) {
  y <- df[units,1,drop=TRUE];
  RawAlleleBFractions(y, x=pos, chromosome=chromosome);
});
print(fracBList);


devSet("fracB");
subplots(length(fracBList), ncol=1);
par(mar=c(3,3,1,1)+0.1);
cols <- as.integer(1+2*callsN[isHet]);
for (kk in seq(along=fracBList)) {
  fracB <- fracBList[[kk]];
  fracBT <- extractSubset(fracB, subset=isHet);
  plot(fracBT, pch=".", col=cols);
  stext(side=3, pos=1, names(fracBList)[kk]);
} # for (kk ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Mirrored and stratified allele B fractions for LOH analysis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
devSet("mFracB");
subplots(length(fracBList), ncol=1);
par(mar=c(3,3,1,1)+0.1);
cols <- as.integer(1+2*callsN[isHet]);
for (kk in seq(along=fracBList)) {
  fracB <- fracBList[[kk]];
  fracBT <- extractSubset(fracB, subset=isHet);
  mFracBT <- extractRawMirroredAlleleBFractions(fracBT);
  plot(mFracBT, pch=".", ylim=c(0,0.6), col=cols);
  stext(side=3, pos=1, names(fracBList)[kk]);
  fit <- segmentByCBS(mFracBT);
  cnr <- extractCopyNumberRegions(fit);
  drawLevels(cnr, col="black");
} # for (kk ...)


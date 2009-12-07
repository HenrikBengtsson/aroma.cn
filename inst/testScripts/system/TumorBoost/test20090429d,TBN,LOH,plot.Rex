if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fullnames translator
fnt <- function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
} # fnt()

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

pattern <- "^(TCGA,OV,testSet,pairs,.*),TBN,.*";
dataSets <- list.files(path=rootPath, pattern=pattern);
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
# Load the normalized tumor data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsTN <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
setFullNamesTranslator(dsTN, fnt);

sampleNames <- sort(unique(getNames(dsTN)));
print(sampleNames);

if (interactive()) {
  sampleName <- textMenu(sampleNames, value=TRUE);
} else {
  sampleName <- sampleNames[1];
}
 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the raw (tumor,normal) data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dataSetOrg <- gsub(pattern, "\\1", dataSet);
ds <- AromaUnitFracBCnBinarySet$byName(dataSetOrg, chipType="*", paths=rootPath);
setFullNamesTranslator(ds, fnt);
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
# Load the genotype call set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Identify available genotype calls
rootPath <- "callData";
rootPath <- Arguments$getReadablePath(rootPath);

pattern <- sprintf("^%s,", getFullName(ds));
dataSets <- list.files(path=rootPath, pattern=pattern);
callTypes <- gsub(pattern, "", dataSets);
# Sanity check
stopifnot(length(callTypes) > 0);

if (interactive()) {
  callType <- textMenu(callTypes, value=TRUE);
} else {
  callType <- callTypes[1];
}
print(callType);

genotypeTag <- callType; 
gsN <- AromaUnitGenotypeCallSet$byName(dataSetOrg, tags=genotypeTag, chipType="*");
setFullNamesTranslator(gsN, fnt);
# Keep only normal genotype files
types <- sapply(gsN, FUN=function(df) getTags(df)[1]);
types <- gsub("[A-Z]$", "", types);
keep <- which(is.element(types, c("10", "11")));
gsN <- extract(gsN, keep);
print(gsN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create an list of matched data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsList <- list(normal=dsN, tumor=dsT, tumorN=dsTN, callsN=gsN);
dsList <- lapply(dsList, FUN=function(ds) {
  idxs <- indexOf(ds, getNames(dsList$normal));
  extract(ds, idxs);
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- getAromaUgpFile(dsList$tumor);
chromosome <- 2;
units <- getUnitsOnChromosome(ugp, chromosome=chromosome);
pos <- getPositions(ugp, units=units);


kk <- 1;

# Extract files of interest
dfList <- lapply(dsList, FUN=getFile, kk);
gf <- dfList[[4]];
dfList <- dfList[-4];

# Stratify on heterozygous normal genotypes
hets <- whichVector(isHeterozygote(gf, units=units, drop=TRUE));

fracBList <- lapply(dfList, FUN=function(df) {
  y <- df[units,1,drop=TRUE];
  RawAlleleBFractions(y, x=pos, chromosome=chromosome);
});
names(fracBList)[3] <- sprintf("%s (normalized using %s calls)", names(fracBList)[3], genotypeTag);
print(fracBList);


devSet(sprintf("fracB,%s", genotypeTag));
subplots(length(fracBList), ncol=1);
par(mar=c(3,3,1,1)+0.1);
for (kk in seq(along=fracBList)) {
  fracB <- fracBList[[kk]];
  fracBT <- extractSubset(fracB, subset=hets);
  plot(fracBT, pch=".");
  stext(side=3, pos=0, sprintf("Heterozygous loci according to %s", genotypeTag));
  stext(side=3, pos=1, names(fracBList)[kk]);
} # for (kk ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Mirrored and stratified allele B fractions for LOH analysis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
devSet(sprintf("mFracB,%s", genotypeTag));
subplots(length(fracBList), ncol=1);
par(mar=c(3,3,1,1)+0.1);
for (kk in seq(along=fracBList)) {
  fracB <- fracBList[[kk]];
  fracBT <- extractSubset(fracB, subset=hets);
  mFracBT <- extractRawMirroredAlleleBFractions(fracBT);
  plot(mFracBT, pch=".", ylim=c(0,0.6));
  stext(side=3, pos=0, sprintf("Heterozygous loci according to %s", genotypeTag));
  stext(side=3, pos=1, names(fracBList)[kk]);
  fit <- segmentByCBS(mFracBT);
  cnr <- extractCopyNumberRegions(fit);
  drawLevels(cnr, col="red", lwd=3);
} # for (kk ...)

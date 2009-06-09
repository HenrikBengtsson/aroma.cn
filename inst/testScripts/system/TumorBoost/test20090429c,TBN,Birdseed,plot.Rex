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
dsTN <- AromaUnitFracBCnBinarySet$byName(dataSet, tags="TBN,Birdseed", chipType="*", paths=rootPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the genotype call set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
genotypeTag <- "Birdseed";
gsN <- AromaUnitGenotypeCallSet$byName(dataSet, tags=genotypeTag, chipType="*");
print(gsN);


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
units <- getUnitsOnChromosome(ugp, 1);
pos <- getPositions(ugp, units=units);


kk <- 1;
dfList <- lapply(dsList, FUN=getFile, kk);
beta <- lapply(dfList, FUN=function(df) df[units,1,drop=TRUE]);
beta <- as.data.frame(beta);
beta <- as.matrix(beta);

subplots(ncol(beta)-1, ncol=1);
par(mar=c(3,3,1,1)+0.1);
ylim <- c(-0.1,1.1);
cols <- as.integer(1+2*beta[,"callsN"]);
for (cc in 1:(ncol(beta)-1)) {
  plot(pos, beta[,cc], pch=".", col=cols, ylim=ylim);
  stext(side=3, pos=1, colnames(beta)[cc]);
}

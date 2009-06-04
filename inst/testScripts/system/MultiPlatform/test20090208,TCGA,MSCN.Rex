if (interactive()) savehistory();
library("aroma.cn");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

tagsList <- list("MSKCC", "Harvard", "Stanford", "Broad");
dataSet <- "TCGA,GBM,testSet,pairs";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load raw CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsList <- lapply(tagsList, FUN=function(tags) {
  AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType="*");
});
# Keep only common samples (just in case)
names <- Reduce(intersect, lapply(dsList, FUN=getNames));
dsList <- lapply(dsList, FUN=extract, names);
print(dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Multi-source CN normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
fitUgp <- AromaUgpFile$byChipType("GenericHuman", tags="100kb");
print(fitUgp);
mscn <- MultiSourceCopyNumberNormalization(dsList, fitUgp=fitUgp);
print(mscn);

dsNList <- process(mscn, verbose=log);
print(dsNList);

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

pattern <- "(^TCGA,OV,testSet,pairs|TumorBoostPaper)";
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
setFullNamesTranslator(ds, fnt);
print(ds);

sampleNames <- sort(unique(getNames(ds)));
print(sampleNames);

if (interactive()) {
  sampleName <- textMenu(sampleNames, value=TRUE);
} else {
  sampleName <- sampleNames[1];
}

flavors <- c("v1", "v2", "v3", "v4");
if (interactive()) {
  flavor <- textMenu(flavors, value=TRUE);
} else {
  flavor <- flavors[4];
}

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
gsN <- AromaUnitGenotypeCallSet$byName(dataSet, tags=genotypeTag, chipType="*");
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
dsList <- list(normal=dsN, tumor=dsT, callsN=gsN);
dsList <- lapply(dsList, FUN=function(ds) {
  idxs <- indexOf(ds, getNames(dsList$normal));
  extract(ds, idxs);
});
print(dsList);

dummy <- lapply(dsList, FUN=function(ds) print(getFile(ds,1)));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Normalize allele B fractions for tumors given matched normals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tbn <- TumorBoostNormalization(dsList$tumor, dsList$normal, gcN=dsList$callsN, flavor=flavor, tags=c("*", genotypeTag));
dsTN <- process(tbn, verbose=log);
print(dsTN);

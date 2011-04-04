if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

pattern <- "(^TCGA,OV,testSet,pairs|TumorBoostPaper)";
pattern <- "^broad.mit.edu_(GBM|OV).Genome_Wide_SNP_6(.*),(.*)";
pattern <- NULL;

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

sampleNames <- sort(unique(getNames(ds)));
print(sampleNames);

if (interactive()) {
  sampleNames <- selectMenu(sampleNames);
} 

flavors <- c("v1", "v2", "v3", "v4")[4];
if (interactive() && length(flavors) > 1) {
  flavor <- textMenu(flavors, value=TRUE);
} else {
  flavor <- flavors[1];
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Load genotype call set
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Identify available genotype calls
rootPath <- "callData";
rootPath <- Arguments$getReadablePath(rootPath);

pattern <- sprintf("^%s,", getFullName(ds));
dataSets <- list.files(path=rootPath, pattern=pattern);
callTypes <- gsub(pattern, "", dataSets);
## Sanity check
stopifnot(length(callTypes) > 0);

if (interactive()) {
  callType <- textMenu(callTypes, value=TRUE);
} else {
  callType <- callTypes[1];
}
print(callType);

genotypeTag <- callType;
gsN <- AromaUnitGenotypeCallSet$byName(dataSet, tags=genotypeTag, chipType="*");
print(gsN);

for (sampleName in sampleNames) {
  pair <- indexOf(ds, sampleName);
  if (length(pair) != 2) {
    verbose && cat(verbose, "Incomplete tumor-normal pair for sample: ", sampleName);
  } else {
    ## Order as (tumor,normal)
    types <- sapply(extract(ds,pair), FUN=function(df) getTags(df)[1]);
    o <- order(types);
    types <- types[o];
    pair <- pair[o];
  
    ## Extract (tumor, normal) pair
    dsPair <- extract(ds, pair);
    dsT <- extract(dsPair, 1);
    dsN <- extract(dsPair, 2);
    
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ## Create a list of matched data sets
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    dsList <- list(normal=dsN, tumor=dsT, callsN=gsN);
    dsList <- lapply(dsList, FUN=function(ds) {
      ds <- setFullNamesTranslator(ds, function(names, ...) {
        pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
        gsub(pattern, "\\1,\\2,\\3", names);
      });
      idxs <- indexOf(ds, getNames(dsList$normal));
      extract(ds, idxs);
    });
    print(dsList);
    
    
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ## Normalize allele B fractions for tumors given matched normals
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    tbn <- TumorBoostNormalization(dsList$tumor, dsList$normal, gcN=dsList$callsN, flavor=flavor, tags=c("*", genotypeTag));
    dsTN <- process(tbn, verbose=log);
    print(dsTN);
  }
}

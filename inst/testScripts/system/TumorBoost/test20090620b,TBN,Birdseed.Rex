if (interactive()) savehistory();
library("aroma.cn");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";

## dataSet <- "TCGA,GBM,onePair";
## sampleName <- "TCGA-12-0620";
dataSet <- "TCGA,OV,testSet,pairs,Broad";
sampleNames <- c("TCGA-04-1335", "TCGA-04-1336", "TCGA-23-1027")

flavor <- c("v1", "v2")[2]

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Load the raw (tumor,normal) data set
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
ds <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
setFullNamesTranslator(ds, function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
});
print(ds);

for (sampleName in sampleNames) {

  pair <- indexOf(ds, sampleName);
  stopifnot(length(pair) == 2);

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
  ## Load the genotype call set
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  genotypeTag <- "Birdseed";
  gsN <- AromaUnitGenotypeCallSet$byName(dataSet, tags=genotypeTag, chipType="*");
  print(gsN);



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
} ## for (sampleName in ...

############################################################################
## HISTORY:
## 2009-06-20
## o Several samples are now processed in a row.
## o Either TBN or TBN,v2 can be used (variable 'flavor').
## o Created from test20090429b,TBN,Birdseed.Rex.
############################################################################

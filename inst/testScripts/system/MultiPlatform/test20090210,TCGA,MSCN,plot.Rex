if (interactive()) savehistory();
library("aroma.cn");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
fig <- 1;

tagsList <- list("MSKCC", "Harvard", "Stanford", "Broad");
tagsList <- lapply(tagsList, FUN=c, "mscn");
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


tags <- Reduce(intersect, lapply(dsList, FUN=getTags));
sites <- sapply(dsList, FUN=function(ds) setdiff(getTags(ds), tags));
platforms <- sapply(dsList, FUN=getPlatform);
chipTypes <- sapply(dsList, FUN=getChipType);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract subset of raw CNs across platforms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Sample of interest
array <- 2;

# Chromosomal region of interest
chromosome <- 3;
region <- c(0,60)*1e6;

cnList <- lapply(dsList, FUN=function(ds) {
  df <- getFile(ds, array);
  extractRawCopyNumbers(df, chromosome=chromosome, region=region);
})
names(cnList) <- sprintf("%s\n%s\n%s", sites, platforms, chipTypes);
print(cnList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Smooth CNs using consecutive bins of width 100kb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xRange <- range(sapply(cnList, FUN=xRange));
cnSList <- lapply(cnList, FUN=function(cn) {
  binnedSmoothing(cn, from=xRange[1], to=xRange[2], by=100e3);
})
print(cnSList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot raw and binned CNs along chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Mlim <- c(-1,1)*3;
xlim <- xRange;

# The name of the tumor/normal pair
name <- getNames(dsList[[1]])[array];

devSet("panels");
layout(seq(along=cnList));
par(mar=c(4.2,4.2,1.3,2.1));
for (kk in seq(along=cnList)) {
  cn <- cnList[[kk]];
  cnS <- cnSList[[kk]];
  plot(cn, xlim=xlim, ylim=Mlim);
  points(cnS, cex=1, col="white");
  points(cnS, cex=0.5, col="red");
  stext(side=3, pos=0, cex=0.8, name);
  stext(side=3, pos=1, cex=0.8, sprintf("Chr%02d (n=%d; s=%d)", 
                               chromosome, nbrOfLoci(cn), nbrOfLoci(cnS)));
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot pairwise binned CNs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("pairs");
Y <- sapply(cnSList, FUN=function(cs) getCNs(cs));
panel <- function(...) { abline(a=0, b=1, lty=2); points(...); }
pairs(Y, pch=20, col="red", lower.panel=NULL, upper.panel=panel, 
                                                     xlim=Mlim, ylim=Mlim);

if (interactive()) savehistory();
library("aroma.cn");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
fig <- 1;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load smoothed CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
rootPath <- "smoothCnData";
rootPath <- Arguments$getReadablePath(rootPath);
dataSets <- list.files(pattern="^TCGA,OV,.*,mscn,v3,.*,by=200kb,.*", path=rootPath);

dsList <- lapply(dataSets, FUN=function(dataSet) {
  AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
});

namesList <- lapply(dsList, FUN=getNames);
uniqueNames <- sort(unique(unlist(namesList, use.names=FALSE)));

# Keep only common samples (just in case)
names <- uniqueNames[1:10];
dsList <- lapply(dsList, FUN=function(ds) {
  idxs <- which(is.element(names, getNames(ds)));
  extract(ds, idxs);
});
dsList <- dsList[(sapply(dsList, length) > 0)];
print(dsList);


tagsList <- lapply(dsList, FUN=getTags);
tags <- Reduce(intersect, tagsList);
sites <- sapply(dsList, FUN=function(ds) setdiff(getTags(ds), tags));
labels <- sapply(sites, FUN=function(ss) paste(ss[1:2], collapse=","));
labels <- gsub(",GenomeWideSNP_6", "(GWS6)", labels);
labels <- gsub(",HG-CGH-244A", "(244A)", labels);
labels <- gsub(",HG-CGH-1x1M", "(1M)", labels);
labels <- gsub(",Human1M-Duo", "(1M)", labels);
sites <- sapply(sites, FUN=function(ss) ss[1]);
platforms <- sapply(dsList, FUN=getPlatform);
chipTypes <- sapply(dsList, FUN=getChipType);
# Sanity check
stopifnot(length(unique(platforms)) == 1);
stopifnot(length(unique(chipTypes)) == 1);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract subset of raw CNs across platforms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Sample of interest
array <- 1;

chromosome <- 6;
region <- NULL;

cnSList <- lapply(dsList, FUN=function(ds) {
  df <- getFile(ds, array);
  extractRawCopyNumbers(df, chromosome=chromosome, region=region);
})
names(cnSList) <- sprintf("%s\n%s\n%s", sites, platforms, chipTypes);
print(cnSList);


xRange <- range(sapply(cnSList, FUN=xRange));

sigmas <- sapply(cnSList, FUN=function(cn) cn$sigma);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot raw and binned CNs along chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Mlim <- c(-1,1)*3;
xlim <- xRange;

# The name of the tumor/normal pair
name <- getNames(dsList[[1]])[array];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Remove shifts
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate global shifts (average over all loci)
yRef <- getSignals(cnSList[[1]]);
deltas <- sapply(cnSList, FUN=function(cn) {
  median(getSignals(cn)-yRef, na.rm=TRUE);
});

cnSM <- sort(Reduce(append, lapply(cnSList, clone)));


siteCols <- c(HMS=1, MSKCC=2, Broad=3, HudsonAlpha=4)
cols <- siteCols[sites];
xScale <- 1e-6;

devSet("smooth;merged;panels");
par(mar=c(2.4,4.2,1.3,2.1));
cn <- cnSM;
cn$sigma <- estimateStandardDeviation(cn);
plot(cn, pch=".", cex=2, col=cols[cn$id], xlim=xlim*xScale, ylim=Mlim, xScale=xScale);
stext(side=3, pos=1/2, cex=0.8, col="red", sprintf("s*=%.3f", cn$sigma));
stext(side=3, pos=0, cex=0.8, name);
stext(side=3, pos=1, cex=0.8, sprintf("Chr%02d (n[s]=%d)", chromosome, nbrOfLoci(cn)));
o <- order(cols);
cols <- cols[o];
labels <- labels[o];
legend("topleft", horiz=TRUE, pch=19, col=cols, labels, bty="n", cex=0.7);
devDone();

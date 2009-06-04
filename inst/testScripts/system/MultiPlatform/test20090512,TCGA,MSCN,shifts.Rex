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

chromosome <- 19;
region <- NULL;

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
  cnS <- binnedSmoothing(cn, from=xRange[1], to=xRange[2], by=100e3);
  cnS$sigma <- estimateStandardDeviation(cnS);
  cnS;
})
print(cnSList);

sigmas <- sapply(cnSList, FUN=function(cn) cn$sigma);
ws <- 1/sigmas^2;
ws <- ws / sum(ws, na.rm=TRUE);

cnList <- mapply(cnList, ws, FUN=function(cn, w) {
  cn$weights <- rep(w, nbrOfLoci(cn));
  list(cn);
})
print(cnList);

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
  plot(cn, pch=".", cex=2, col="gray", xlim=xlim, ylim=Mlim);
  points(cnS, cex=1, col="white");
  points(cnS, cex=0.5, col="red");
  stext(side=3, pos=1/2, cex=0.8, col="red", sprintf("s*=%.3f => w=%.2f", cnS$sigma, cn$weights[1]));
  stext(side=3, pos=0, cex=0.8, name);
  stext(side=3, pos=1, cex=0.8, sprintf("Chr%02d (n=%d; n[s]=%d)", 
                               chromosome, nbrOfLoci(cn), nbrOfLoci(cnS)));
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Remove shifts
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate global shifts (average over all loci)
yRef <- getSignals(cnSList[[1]]);
deltas <- sapply(cnSList, FUN=function(cn) {
  median(getSignals(cn)-yRef, na.rm=TRUE);
});

# Subtract shifts from signals
cnTList <- mapply(cnList, deltas, FUN=function(cn, delta) {
  cn <- clone(cn);
  cn$y <- cn$y - delta;
  list(cn);
});

cnTSList <- lapply(cnTList, FUN=function(cn) {
  cnS <- binnedSmoothing(cn, from=xRange[1], to=xRange[2], by=100e3);
  cnS$sigma <- estimateStandardDeviation(cnS);
  cnS;
})

cnSM <- sort(Reduce(append, lapply(cnSList, clone)));
cnTSM <- sort(Reduce(append, lapply(cnTSList, clone)));

devSet("smooth;merged;panels");
layout(1:2);
par(mar=c(4.2,4.2,1.3,2.1));
cn <- cnSM;
cn$sigma <- estimateStandardDeviation(cn);
plot(cn, pch=".", cex=2, col=cn$id, xlim=xlim, ylim=Mlim);
stext(side=3, pos=1/2, cex=0.8, col="red", sprintf("s*=%.3f", cn$sigma));
stext(side=3, pos=0, cex=0.8, name);
stext(side=3, pos=1, cex=0.8, sprintf("Chr%02d (n[s]=%d)", chromosome, nbrOfLoci(cn)));

cn <- cnTSM;
cn$sigma <- estimateStandardDeviation(cn);
plot(cn, pch=".", cex=2, col=cn$id, xlim=xlim, ylim=Mlim);
stext(side=3, pos=1/2, cex=0.8, col="red", sprintf("s*=%.3f (shifted)", cn$sigma));
stext(side=3, pos=0, cex=0.8, name);
stext(side=3, pos=1, cex=0.8, sprintf("Chr%02d (n[s]=%d)", chromosome, nbrOfLoci(cn)));

devDone();

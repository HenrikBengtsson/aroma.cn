if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

rootPath <- "rawCnData";
rootPath <- Arguments$getReadablePath(rootPath);
dataSets <- list.files(path=rootPath);

# Sanity check
stopifnot(length(dataSets) > 0);

if (interactive()) {
  dataSet <- textMenu(dataSets, value=TRUE);
} else {
  dataSet <- dataSets[1];
}
print(dataSet);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup CN data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*");
print(ds);

ugp <- getAromaUgpFile(ds);
print(ugp);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Bin to 200kb uniformly distributed target loci
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
targetUgp <- AromaUgpFile$byChipType("GenericHuman", tags="200kb");
print(targetUgp);

sm <- TotalCnBinnedSmoothing(ds, targetUgp=targetUgp);
print(sm);

dsS <- process(sm, verbose=verbose);
print(dsS);


# Plot an example
ii <- 1;
df <- getFile(ds, ii);
dfS <- getFile(dsS, ii);
cn <- extractRawCopyNumbers(df, logBase=2, chromosome=4);
cnS <- extractRawCopyNumbers(dfS, logBase=NULL, chromosome=4);
toPNG(getName(dfS), tags=getTags(dsS), width=840, aspectRatio=0.6, {
  plot(cn, ylim=c(-5,+5), col="#666666", pch=".", xlim=c(105,120));
  points(cnS, col="#ff0000");
  stext(side=3, pos=0, getName(dfS));
  stext(side=4, pos=0, cex=0.8, sprintf("Chip type: %s", getChipType(df)));
  stext(side=4, pos=1, cex=0.8, sprintf("Target chip type: %s", getChipType(dfS)));
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Bin to target chip type Mapping10K_Hind142
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
targetUgp <- AromaUgpFile$byChipType("Mapping10K_Xba142");
print(targetUgp);

sm <- TotalCnBinnedSmoothing(ds, targetUgp=targetUgp, tags=c("*", "10K"));
print(sm);

dsS <- process(sm, verbose=verbose);
print(dsS);


# Plot an example
ii <- 1;
df <- getFile(ds, ii);
dfS <- getFile(dsS, ii);
cn <- extractRawCopyNumbers(df, logBase=2, chromosome=4);
cnS <- extractRawCopyNumbers(dfS, logBase=NULL, chromosome=4);
toPNG(getName(dfS), tags=getTags(dsS), width=840, aspectRatio=0.6, {
  plot(cn, ylim=c(-5,+5), col="#666666", pch=".", xlim=c(105,120));
  points(cnS, col="#ff0000");
  stext(side=3, pos=0, getName(dfS));
  stext(side=4, pos=0, cex=0.8, sprintf("Chip type: %s", getChipType(df)));
  stext(side=4, pos=1, cex=0.8, sprintf("Target chip type: %s", getChipType(dfS)));
});

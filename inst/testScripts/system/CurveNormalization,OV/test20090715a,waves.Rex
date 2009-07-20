if (interactive()) savehistory();
library("aroma.cn");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tcgaFNT01 <- function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scan for data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

dataSet <- "TCGA,OV,waves,HudsonAlpha";
dataSet <- "TCGA,OV,waves,Broad";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsR <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
setFullNamesTranslator(dsR, tcgaFNT01);
print(dsR);

ugp <- getAromaUgpFile(dsR);
print(ugp);

sampleNames <- sort(unique(getNames(dsR)));
idxs <- lapply(sampleNames, FUN=function(name) {
  idxs <- indexOf(dsR, name);
  o <- order(getFullNames(dsR)[idxs]);
  idxs[o];
});
ns <- sapply(idxs, FUN=length);
stopifnot(all(ns == 2));
idxs <- matrix(unlist(idxs, use.names=FALSE), nrow=2, byrow=FALSE);
rownames(idxs) <- c("T", "N");
colnames(idxs) <- sampleNames;
print(idxs);

setFullNamesTranslator(dsR, NULL);
ds <- extract(dsR, idxs["T",]);
dsTarget <- extract(dsR, idxs["N",]);
pcn <- XYCurveNormalization(ds, dsTarget);
pcn <- PrincipalCurveNormalization(ds, dsTarget);
print(pcn);

dsN <- process(pcn, verbose=verbose);
print(dsN);

chr <- 7;
units <- getUnitsOnChromosome(ugp, chromosome=chr);
pos <- getPositions(ugp, units=units);

lim <- c(0,6);
Alim <- c(-5,3);
Mlim <- c(-3,3);

fullnames <- getFullNames(dsR);
fullnames <- fullnames[idxs[,1]];
dsList <- list(dsR, dsN);
subplots(3*length(dsList), ncol=3, byrow=TRUE);
for (ds in dsList) {
  dsPair <- extract(ds, indexOf(ds, fullnames));
  theta <- extractMatrix(dsPair, units=units);

  plot(theta[,2:1], pch=".", xlim=lim, ylim=lim);
  abline(a=0, b=1, lwd=2, col="#999999");
  subset <- sample(1:nrow(theta), size=nrow(theta)/8);
  fit <- fitPrincipalCurve(theta[subset,]);
  xy <- fit$s[,2:1];
  xy <- xy[order(xy[,1]),];
  lines(xy, col="red", lwd=2);

  plotMvsA(theta, pch=".", Alim=Alim);
  abline(h=0, lwd=2, col="#999999");

  cn <- RawCopyNumbers(cn=log2(theta[,2]/theta[,1]), x=pos);
  plot(cn, pch=".");
  stext(side=3, pos=1, sprintf("Chr %d", chr));
  cnS <- binnedSmoothing(cn, by=500e3);
  points(cnS, pch=".", cex=3, col="red");
} # for (ds ...)


############################################################################
# HISTORY:
# 2009-07-15
# o Created.
############################################################################

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

figPath <- Arguments$getWritablePath("figures");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

pattern <- "^(TCGA,OV,testSet,pairs,[^,]*),(.*),TBN(,|)(v[0-9]|),.*";
dataSets <- list.files(path=rootPath, pattern=pattern);
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

flavor <- gsub(pattern, "\\3", dataSet);
if (flavor == "") {
  flavor <- "v1";
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the normalized tumor data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsTN <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
setFullNamesTranslator(dsTN, fnt);

sampleNames <- sort(unique(getNames(dsTN)));
print(sampleNames);

if (interactive()) {
  sampleName <- textMenu(sampleNames, value=TRUE);
} else {
  sampleName <- sampleNames[1];
}
 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the raw (tumor,normal) data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dataSetOrg <- gsub(pattern, "\\1,\\2", dataSet);
ds <- AromaUnitFracBCnBinarySet$byName(dataSetOrg, chipType="*", paths=rootPath);
setFullNamesTranslator(ds, fnt);
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


siteTag <- getTags(ds);
siteTag <- siteTag[4];



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the genotype call set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Identify available genotype calls
rootPath <- "callData";
rootPath <- Arguments$getReadablePath(rootPath);

pattern <- sprintf("^%s,", dataSetOrg);
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
gsN <- AromaUnitGenotypeCallSet$byName(dataSetOrg, tags=genotypeTag, chipType="*");
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
dsList <- list(normal=dsN, tumor=dsT, tumorN=dsTN, callsN=gsN);
rm(dsN, dsT, dsTN, gsN);
dsList <- lapply(dsList, FUN=function(ds) {
  idxs <- indexOf(ds, getNames(dsList$normal));
  extract(ds, idxs);
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the raw CN data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);
ds <- AromaUnitTotalCnBinarySet$byName(dataSetOrg, chipType="*", paths=rootPath);
setFullNamesTranslator(ds, fnt);
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- getAromaUgpFile(dsList$tumor);

if (interactive()) {
  chromosome <- textMenu(getChromosomes(ugp), value=TRUE);
} else {
  chromosome <- 2;
}
chrTag <- sprintf("Chr%02d", chromosome);
units <- getUnitsOnChromosome(ugp, chromosome=chromosome);

# Identify SNPs only
platform <- getPlatform(ugp);
if (platform == "Affymetrix") {
  require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix");
  snpPattern <- "^SNP";
} else if (platform == "Illumina") {
  snpPattern <- "^rs[0-9]";
} else {
  throw("Unknown platform: ", platform);
}
unf <- getUnitNamesFile(ugp);
unitNames <- getUnitNames(unf, units=units);

# Identify SNP units
keep <- (regexpr(snpPattern, unitNames) != -1);
units <- units[keep];

pos <- getPositions(ugp, units=units);

# Extract total CNs
C <- extractMatrix(dsPair, units=units);
C <- 2*C[,1]/C[,2];

# Extract Allele B fractions
kk <- 1;
dfList <- lapply(dsList, FUN=getFile, kk);
beta <- lapply(dfList, FUN=function(df) df[units,1,drop=TRUE]);
beta <- as.data.frame(beta);
beta <- as.matrix(beta);
names <- colnames(beta);
names[names == "tumorN"] <- "normalized tumor";

# Plot dimensions
x <- pos/1e6;
xlim <- range(x, na.rm=TRUE);
xlab <- "Position (Mb)";

width <- 840;
width <- 1280;
aspect <- 0.6*1/3;

# Plot total CNs
ylim <- c(0,6);
ylab <- "Copy number";
figName <- sprintf("%s,%s,%s,CN", siteTag, sampleName, chrTag);
pathname <- filePath(figPath, sprintf("%s.png", figName));
if (!isFile(pathname)) {
  fig <- devNew("png", pathname, label=figName, width=width, height=aspect*width);
  par(mar=c(2.7,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
  par(mar=c(1.7,2.5,1.1,1)+0.1);
  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=c(0,2,4,6));
  points(x, C, pch=".");
  label <- sprintf("%s", sampleName);
  stext(side=3, pos=0, label);
  stext(side=3, pos=1, chrTag);
  devDone();
}

# Plot Allele B fractions
ylim <- c(-0.05,1.05);
ylim <- c(-0.1,1.1);
ylab <- "Allele B Fraction";
cols <- as.integer(beta[,"callsN"] != 1) + 1L;

for (cc in 1:3) {
  tag <- colnames(beta)[cc];
  name <- names[cc];
  figName <- sprintf("%s,%s,%s,%s,fracB", siteTag, sampleName, tag, chrTag);
  pathname <- filePath(figPath, sprintf("%s.png", figName));
  if (!isFile(pathname)) {
    fig <- devNew("png", pathname, label=figName, width=width, height=aspect*width);
    par(mar=c(2.7,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
    par(mar=c(1.7,2.5,1.1,1)+0.1);
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE);
    axis(side=1);
    axis(side=2, at=c(0,1/2,1));
    points(x, beta[,cc], pch=".", col=cols);
    label <- sprintf("%s (%s)", sampleName, name);
    stext(side=3, pos=0, label);
    stext(side=3, pos=1, chrTag);
    devDone();
  }
}

for (tag in c("tumor", "tumorN")) {
  figName <- sprintf("%s,%s,%s,%s,betaVsBeta", siteTag, sampleName, tag, chrTag);
  pathname <- filePath(figPath, sprintf("%s.png", figName));
  if (!isFile(pathname)) {
    fig <- devNew("png", pathname, label=figName, width=640, height=640);
    par(mar=c(2.7,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
    lim <- c(-0.1,1.1);
    plot(beta[,c("normal", tag)], pch=".", xlim=lim, ylim=lim);
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, chrTag);
    devDone();
  }
}


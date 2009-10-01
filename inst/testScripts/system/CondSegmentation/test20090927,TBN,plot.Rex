if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");

kappa <- 0.6;
#kappa <- 0.98;

figPath <- Arguments$getWritablePath("figures");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

pattern <- "^(TCGA,OV,testSet,pairs,[^,]*),(.*),TBN(,|)(v[0-9]|),.*";
#pattern <- "^(TCGA,OV,testSet,pairs,[^,]*),(.*)(,TBN,.*|)$";
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the normalized tumor data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsTN <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*");
setFullNamesTranslator(dsTN, function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
}); 

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
ds <- AromaUnitFracBCnBinarySet$byName(dataSetOrg, chipType="*");
setFullNamesTranslator(ds, function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
});
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
print(gsN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create an list of matched data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsList <- list(normal=dsN, tumor=dsT, tumorN=dsTN, callsN=gsN);
rm(dsN, dsT, dsTN, gsN);
dsList <- lapply(dsList, FUN=function(ds) {
  ds <- setFullNamesTranslator(ds, function(names, ...) {
    pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
    gsub(pattern, "\\1,\\2,\\3", names);
  });
  idxs <- indexOf(ds, getNames(dsList$normal));
  extract(ds, idxs);
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load the raw CN data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ds <- AromaUnitTotalCnBinarySet$byName(dataSetOrg, chipType="*");
setFullNamesTranslator(ds, function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
});
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
# Segment the total CNs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Clim <- c(0,5);
rhoLim <- c(0,1);

C1lab <- expression(C[minor]);
C2lab <- expression(C[major]);
Clab <- expression(C);
rholab <- expression(rho);
  

dataList <- list(
  dfCallN = getFile(dsList$callsN, 1),
  dfThetaT = getFile(dsPair, 1),
  dfThetaN = getFile(dsPair, 2),
  dfBetaT = getFile(dsList$tumor, 1),
  dfBetaTN = getFile(dsList$tumorN, 1)
);

chromosomes <- 1:24;
sampleName <- getName(dataList[[1]]);
chipType <- getChipType(dataList[[1]]);

key <- list(dataSetNames=sapply(dsList, FUN=getFullName), chipType=chipType, sampleName=sampleName, chromosomes=chromosomes);
dirs <- c("aroma.cn", "TumorPurifier");

data <- loadCache(key=key, dirs=dirs);
if (is.null(data)) {
  verbose && enter(verbose, "Identifying all heterzygous SNPs");
  isHet <- isHeterozygous(dataList$dfCallN);
  hets <- whichVector(isHet);
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, hets);
  verbose && exit(verbose);
  
  figName <- sprintf("C1C2,kappa=%.2f", kappa);
  plot(NA, pch=".", xlim=Clim, ylim=Clim, xlab=C1lab, ylab=C2lab);
  abline(a=0, b=1, lty=2);
  stext(side=3, pos=0, sampleName);
  
  data <- NULL;
  for (chr in chromosomes) {
    verbose && enter(verbose, sprintf("Chromosome %d", chr));
  
    verbose && enter(verbose, "Extracting total copy numbers");
    thetaT <- extractRawGenomicSignals(dataList$dfThetaT, chromosome=chr);
    thetaN <- extractRawGenomicSignals(dataList$dfThetaN, chromosome=chr);
    cn <- divideBy(thetaT, thetaN);
    cn$y <- 2*cn$y;
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Segmenting");
    fit <- segmentByCBS(cn, cache=TRUE);
    verbose && print(verbose, fit);
    regions <- extractCopyNumberRegions(fit);
    verbose && print(verbose, regions);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Extracting tumor allele B fractions");
    fracB <- extractRawGenomicSignals(dataList$dfBetaTN, chromosome=chr, units=hets);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Calculating average |beta| for heterozygotes region by region");
    dataCC <- as.data.frame(regions);
    for (rr in seq(length=nrow(dataCC))) {
      reg <- c(dataCC[rr,"start"], dataCC[rr,"stop"]);
      fracBRR <- extractRegion(fracB, region=reg);
      n <- nbrOfLoci(fracBRR);
      if (n > 0) {
        rho <- abs(fracBRR$y-1/2);
        avg <- mean(rho, na.rm=TRUE);
      } else {
        avg <- as.double(NA);
      }
      dataCC[rr,"meanRho"] <- avg;
      dataCC[rr,"nRho"] <- n;
    } # for (rr ...)
    verbose && exit(verbose);
  
    C <- dataCC[,"mean"];
    rho <- dataCC[,"meanRho"];
    C1 <- (1/2-rho)*C;
    C2 <- (1/2+rho)*C;
    dataCC[,"C1"] <- C1;
    dataCC[,"C2"] <- C2;
    data <- rbind(data, dataCC);
  
    xy <- cbind(x=C1,y=C2);
    points(xy, pch=".");
    text(xy[,1],xy[,2], dataCC$chromosome);
    verbose && exit(verbose);
  } # for (chr ...)
  saveCache(data, key=key, dirs=dirs);
  } # if (is.null(data))

C1t <- c(0, 0, 1, 0, 1);
C2t <- c(0, 1, 1, 2, 2);
C1x <- kappa*C1t + (1-kappa)*1;
C2x <- kappa*C2t + (1-kappa)*1;

nCount <- 150;
dataT <- subset(data, count > nCount);
dataT <- data;
C <- dataT[,"mean"];
rho <- dataT[,"meanRho"];
C1 <- (1/2-rho)*C;
C2 <- (1/2+rho)*C;

cex <- log10(data[,"count"]/10);

figName <- sprintf("C1C2,nCount=%.0f", nCount);
devSet(figName);
plot(NA, xlim=Clim, ylim=Clim, xlab=C1lab, ylab=C2lab);
points(C1, C2, pch=".", cex=cex);
stext(side=3, pos=0, sampleName);
stext(side=4, pos=0, cex=0.7, chipType);
stext(side=4, pos=1, cex=0.7, sprintf("At least %.0f loci in TCN", nCount));
text(C1,C2, labels=data$chromosome, cex=0.8);
abline(a=0, b=1, lty=2);

figName <- sprintf("C1C2,nCount=%.0f,kappa=%.2f", nCount, kappa);
devSet(figName);
plot(NA, xlim=Clim, ylim=Clim, xlab=C1lab, ylab=C2lab);
points(C1x,C2x, pch=19, cex=2, col="red");
points(C1, C2, pch=".", cex=cex);
stext(side=3, pos=0, sampleName);
stext(side=4, pos=0, cex=0.7, chipType);
stext(side=3, pos=1, sprintf("kappa=%.2f", kappa));
stext(side=4, pos=1, cex=0.7, sprintf("At least %.0f loci in TCN", nCount));
text(C1,C2, labels=data$chromosome, cex=0.8);
abline(a=0, b=1, lty=2);


stop();

figName <- sprintf("Crho");
devSet(figName);
plot(NA, xlim=Clim, ylim=rhoLim, xlab=expression(C), ylab=expression(rho[T]));
stext(side=3, pos=0, sampleName);
stext(side=4, pos=0, cex=0.7, chipType);
stext(side=4, pos=1, cex=0.7, "At least 200 loci in TCN");
points(xy, pch=".");
text(xy, labels=data$chromosome);


if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");

purify <- FALSE;

kappa <- 0.5;
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

getNbrOfModes <- function(y, choices=1:2, default=0) {
  nbrOfModes <- default;
  if ((length(y)>1) && require(mclust)) {
    emFit <- Mclust(y, G=choices, model="E");
    nbrOfModes <- emFit$G;
    str(emFit);
  }
  nbrOfModes;
}

estimatePloidy <- function(C, weights=NULL, na.rm=TRUE) {
  1/weightedMedian(C, weights, na.rm=na.rm)
}
  
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
Clim <- c(0,4);
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
dirs <- c("aroma.cn", "TumorPurifier", "test");

data <- loadCache(key=key, dirs=dirs);
## data <- NULL;
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
    cn$y <- log2(cn$y);
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
    
    verbose && enter(verbose, "Extracting total copy numbers for hets only");
    thetaT <- extractRawGenomicSignals(dataList$dfThetaT, chromosome=chr, units=hets);
    thetaN <- extractRawGenomicSignals(dataList$dfThetaN, chromosome=chr, units=hets);
    cnH <- divideBy(thetaT, thetaN);
    cnH$y <- 2*cnH$y;
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating average major and minor CNs for hets region by region");
    dataCC <- as.data.frame(regions);
    dataCC$mean <- 2*2^dataCC$mean;
    
    for (rr in seq(length=nrow(dataCC))) {
      reg <- c(dataCC[rr,"start"], dataCC[rr,"stop"]);
      fracBRR <- extractRegion(fracB, region=reg);
      cnHRR <- extractRegion(cnH, region=reg);
      n <- nbrOfLoci(fracBRR);
      stopifnot(n==nbrOfLoci(cnHRR));
      AI <- NA
      if (n > 0) {
        rho <- fracBRR$y-1/2;
        nbrOfModes <- getNbrOfModes(fracBRR$y);
        AI <- TRUE;
        if (nbrOfModes==1) {
          AI <- FALSE;
          cat(verbose, "Only one mode in rho for region ", reg);
        } else {
          rho <- abs(rho);
        }
        C1 <- (1/2-rho)*cnHRR$y;
        C2 <- (1/2+rho)*cnHRR$y;
        avg <- mean(rho, na.rm=TRUE);
        avgC1 <- mean(C1, na.rm=TRUE);
        avgC2 <- mean(C2, na.rm=TRUE);
      } else {
        avg <- as.double(NA);
        avgC1 <- as.double(NA);
        avgC2 <- as.double(NA);
      }
      dataCC[rr,"meanRho"] <- avg;
      dataCC[rr,"nRho"] <- n;
      dataCC[rr,"C1"] <- avgC1;
      dataCC[rr,"C2"] <- avgC2;
      dataCC[rr,"AI"] <- AI;
    } # for (rr ...)
    verbose && exit(verbose);
  
    data <- rbind(data, dataCC);
  
    xy <- cbind(x=dataCC[,"C1"],y=dataCC[,"C2"]);
    points(xy, pch=".");
    text(xy[,1],xy[,2], dataCC$chromosome);
    verbose && exit(verbose);
  } # for (chr ...)
  saveCache(data, key=key, dirs=dirs);
  } # if (is.null(data))

C1t <- c(0, 0, 1, 0, 1);
C2t <- c(0, 1, 1, 2, 2);
if (!purify) {
  C1t <- kappa*C1t + (1-kappa)*1;
  C2t <- kappa*C2t + (1-kappa)*1;
}

## estimate ploidy :)
verbose && enter(verbose, "Estimating ploidy")
isAB <- whichVector(!data$AI); ## "allelic balance"
ploidy <- 2*2/weightedMedian(data$mean[isAB], data$nRho[isAB], na.rm=TRUE);
## NB 2*2 is weird but true I think
verbose && cat(verbose, "ploidy: ", ploidy)
verbose && exit(verbose)

verbose && enter(verbose, "Correcting copy numbers from ploidy")
toCorrect <- c("mean", "C1", "C2");
for (tt in toCorrect) {
  data[, tt] <- data[, tt]*ploidy/2;
}
verbose && exit(verbose)

if (purify) {
  verbose && enter(verbose, "Correcting copy numbers from purity")
  data[, "mean"] <- (data[, "mean"]-2*(1-kappa))/kappa;
  data[, "C1"] <- (data[, "C1"]-(1-kappa))/kappa;
  data[, "C2"] <- (data[, "C2"]-(1-kappa))/kappa;
  verbose && exit(verbose)
}

nCount <- 50;
dataT <- subset(data, nRho > nCount);
## dataT <- data;
C <- dataT[,"mean"];
rho <- dataT[,"meanRho"];
C1 <- dataT[,"C1"];
C2 <- dataT[,"C2"];
## isAB <- whichVector(!data$AI); ## "allelic balance"
## C1[isAB] <- C[isAB]/2;
## C2[isAB] <- C[isAB]/2;

enter(verbose, "Defining weights for geneplotter")
dataR <- dataT[, c("nRho", "C1", "C2")];
weight <- dataT$nRho%/%nCount;
C1GP <- rep(dataT$C1, times=weight);
C2GP <- rep(dataT$C2, times=weight);
exit(verbose);

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
points(C1t,C2t, pch=19, cex=2, col="red");
points(C1, C2, pch=".", cex=cex);
stext(side=3, pos=0, sampleName);
stext(side=4, pos=0, cex=0.7, chipType);
stext(side=3, pos=1, sprintf("kappa=%.2f", kappa));
stext(side=4, pos=1, cex=0.7, sprintf("At least %.0f loci in TCN", nCount));
text(C1,C2, labels=data$chromosome, cex=0.8);
abline(a=0, b=1, lty=2);

cex <- 2;
if (require(geneplotter)) {
  figName <- sprintf("C1C2,nCount=%.0f,smoothScatter", nCount);
  devSet(figName);
  Clim <- c(0, ceiling(max(C2)));
  fileName <- sprintf("%s,%s.png", figName, sampleName);
  pathName <- file.path(figPath, fileName);
  png(pathName, width=600, height=600);
  par(mar=c(2.5,2.4,1.3,1)+0.1, mgp=c(1.4,0.3,0), cex=cex);
  smoothScatter(C1GP, C2GP, xlim=Clim, ylim=Clim, xlab=C1lab, ylab=C2lab);
##   stext(side=3, pos=0, sampleName);
  stext(side=4, pos=0, sprintf("%s,%s", sampleName, chipType), cex=1);
  stext(side=3, pos=1, sprintf("CNs corrected for ploidy (estimated as %.2f)", ploidy));

  stext(side=4, pos=1, sprintf(">%.0f heterozygous loci in each region", nCount), cex=1);
  abline(a=0, b=1, lty=2);
  if (purify) {
    abline(h=seq(Clim[1], Clim[2]), lty=4, col="#999999");
    abline(v=seq(Clim[1], Clim[2]), lty=4, col="#999999");
  }
  devDone();
}

############################################################################
## HISTORY:
## 2009-10-14
## o Added an estimate of ploidy and correction for ploidy.
## o Added function 'getNbrOfModes' to test unimodality of BAF.
## When one mode only (no AI), rho is estimated as mean(BAF)-1/2
## (ie essentially 0) instead of mean(abs(BAF-1/2)).
## o Created from test20090927.Rex.
############################################################################

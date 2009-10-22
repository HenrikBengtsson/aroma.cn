if (interactive()) savehistory();
library("aroma.affymetrix");
library("aroma.cn");
library("R.menu");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
extractTcgaPairs <- function(ds, ...) {
  library("aroma.tcga");
  library("gsubfn");
  ps <- BiospecimenCoreResource$getBarcodePatterns();
  pattern <- ps$aliqoutBarcode;
  setFullNamesTranslator(ds, function(names, ...) {
    strapply(names, pattern=pattern, FUN=function(...) {
      x <- unlist(list(...)); n <- length(x);
      sprintf("%s,%s", x[2], x[6]);
    });
  });

  # tumor-normal pairs
  names <- getNames(ds);
  t <- table(names);
  keep <- names[t == 2];
  dsT <- extract(ds, unique(keep));
  names <- getFullNames(dsT);
  dsT <- extract(dsT, order(names));

  dsT;
} # extractTcgaPairs()


figPath <- Arguments$getWritablePath("figures");

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
fig <- 1;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
rootPath <- Arguments$getReadablePath("totalAndFracBData");
dataSets <- list.files(path=rootPath);
dataSet <- textMenu(dataSets, value=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load total DNA signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsT <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*");
print(dsT);

# Assert that annotation data exists
ugp <- getAromaUgpFile(dsT);
print(ugp);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load allele B fractions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsB <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*");
print(dsB);

dsList <- list(total=dsT, fracB=dsB);

if (getName(dsList$total) == "TCGA") {
  dsList <- lapply(dsList, FUN=extractTcgaPairs);
}

print(dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select sample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ds <- dsList$total;
sampleNames <- unique(getNames(ds));
sampleName <- textMenu(sampleNames, value=TRUE);

dsListT <- lapply(dsList, FUN=function(ds) {
  extract(ds, sampleName);
});
print(dsListT);

extractTcgaSignals <- function(dsList, chromosome, ...) {
  ugp <- getAromaUgpFile(dsList$total);
  units <- getUnitsOnChromosome(ugp, chromosome=chromosome);
  ds <- dsList$total;
  dfT <- getFile(ds, 1);
  dfN <- getFile(ds, 1);
  theta <- extractRawCopyNumbers(df, chromosome=chromosome, units=units);
} # extractTcgaSignals()

dfList <- lapply(dsListT, FUN=function(ds) {
  df <- getFile(ds, indexOf(ds, sampleName));
  print(df);

# Is data log2-ed?
isLog2 <- hasTag(df, "log2ratio");

tags <- getTags(df);
tags <- setdiff(tags, c("log2ratio", "total"));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chromosomes <- getChromosomes(ugp);
chromosome <- textMenu(chromosomes, value=TRUE);
chrTag <- sprintf("Chr%02d", chromosome);

region <- textRegion();
if (is.null(region)) {
  regionTag <- NULL;
} else {
  region <- 1e6 * region;
  regionTag <- sprintf("%.2f-%.2fMb", region[1]/1e6, region[2]/1e6);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cn <- extractRawCopyNumbers(df, chromosome=chromosome, region=region);
if (isLog2) {
  cn$y <- 2*2^cn$y;
}
print(cn);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select region interactively
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
pa <- textPlotAnnotations(pa="default");
attachLocally(pa);
cnT <- clone(cn);
cnS <- binnedSmoothing(cnT, by=200e3);
replot <- TRUE;
doPNG <- FALSE;

while (TRUE) {
  print(cnT);
  paTag <- digest(pa, algo="crc32");
  figName <- paste(c(sampleName, chrTag, regionTag, paTag), collapse=",");
  cat(figName);

  if (replot) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Plot total CN track
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (doPNG) {
      filename <- sprintf("%s,TCN.png", figName);
      pathname <- file.path(figPath, filename);
      devNew("png", pathname, width=width, height=height);
    } else {
      devSet(figName);
    }
    par(mar=c(3.5,3.5,1.3,1.1), mgp=c(1.8, 0.8, 0));
    par(cex.axis=1.5, cex.lab=1.5);
    plot(cnT, col=rawCol, xlab=xlab, ylab=ylab, ylim=Clim, xScale=xScale);
    if (smoothCol != "NA") {
      points(cnS, cex=1.5*smoothCex, col="white", xScale=xScale);
      points(cnS, cex=smoothCex, col=smoothCol, xScale=xScale);
    }
    if (showSampleName) {
      stext(side=3, pos=0, sampleName, cex=1.5);
    }
    if (showChr) {
      stext(side=3, pos=1, chrTag, cex=1.5);
    }
    devDone();
    replot <- FALSE;
  }

  choices <- c("Change chromosome", "Change region", "Settings", "Write PNG");
  ans <- textMenu(choices);
  if (ans == 1) {
    chromosome0 <- chromosome;
    chromosome <- textMenu(chromosomes, value=TRUE);
    chrTag <- sprintf("Chr%02d", chromosome);
    cn <- extractRawCopyNumbers(df, chromosome=chromosome);
    if (isLog2) {
      cn$y <- 2*2^cn$y;
    }
    cnT <- clone(cn);
    cnS <- binnedSmoothing(cnT, by=200e3);
    region <- NULL;
    regionTag <- NULL;
    replot <- TRUE;
  } else if (ans == 2) {
    if (!is.null(region)) {
      region <- region / 1e6;
    }
    region0 <- region;
    region <- textRegion(region=region);
    region <- 1e6 * region;
    if (!isTRUE(all.equal(region, region0))) {
      regionTag <- sprintf("%.2f-%.2fMb", region[1]/1e6, region[2]/1e6);
      cnT <- extractRegion(cn, region=region);
      cnS <- binnedSmoothing(cnT, by=200e3);
      replot <- TRUE;
    }
  } else if (ans == 3) {
    pa0 <- pa;
    pa <- textPlotAnnotations(pa);
    if (!isTRUE(all.equal(pa, pa0))) {
      attachLocally(pa);
      str(pa);
      replot <- TRUE;
    }
  } else if (ans == 4) {
    replot <- TRUE;
  }
} # while()


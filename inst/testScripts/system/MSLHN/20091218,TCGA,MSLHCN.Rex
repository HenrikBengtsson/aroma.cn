if (interactive()) savehistory();
library("aroma.cn");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
extractMergedLHTrack <- function(lhList, chromosome, sort=TRUE, ...) {
  chrTag <- sprintf("chr%02d", chromosome);
  lhs <- lapply(lhList, FUN=function(lhs) lhs[[chrTag]]);

  # Merge the all LH tracks
  lh <- clone(lhs[[1]]);
  lhs <- lhs[-1];
  for (kk in seq(along=lhs)) {
    lh <- append(lh, lhs[[kk]]);
  }

  if (sort) {
    lh <- sort(lh);
  }

  lh;
} # extractMergedLHTrack()

plotLHTracks <- function(lhList, ..., lhLim=c(0,1), g=NULL, h=NULL) {
  subplots(seq(along=chromosomes), ncol=3);
  par(mar=c(2,2,0,0.5)+0.1, mgp=c(3,0.7,0));

  chrTags <- names(lhList[[1]]);
  chromosomes <- gsub("chr", "", chrTags);
  chromosomes <- as.integer(chromosomes);
  for (chr in chromosomes) {
    chrTag <- sprintf("chr%02d", chr);
    lh <- extractMergedLHTrack(lhList, chromosome=chr);
    if (nbrOfLoci(lh) < 2) next;
    # Transform signals
    if (!is.null(g)) {
      lh$y <- g(lh$y);
    }
    plot(lh, col=lh$id, ylim=lhLim);
    stext(side=3, pos=1, line=-1, chrTag);
    stext(side=3, pos=0, line=-1, sprintf("col=%d", chr), col=chr);

    if (is.numeric(h)) {
      abline(h=h, lty=3, col="#999999");
    }

  } # for (chr ...)
} # plotLHTracks()

stop()

figPath <- Arguments$getWritablePath("figures/MSLHN/");

tcgaSampleFnt <- function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
}


log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "TCGA,OV";

tbn <- FALSE;

if (tbn) {
  figPath <- file.path(figPath, "TBN");
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load BAF data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tagsList <- list("BeadStudio,XY", "CRMAv2");
dsList <- lapply(tagsList, FUN=function(tags) {
  # Get TumorBoost normalized data?
  if (tbn) {
    tags <- c(tags, "TBN,NGC");
  }
  ds <- AromaUnitFracBCnBinarySet$byName(dataSet, tags=tags, chipType="*");
  setFullNamesTranslator(ds, tcgaSampleFnt);

  # Keep only tumors
  types <- sapply(ds, FUN=function(df) getTags(df)[1]);
  keep <- grep("^01", types);
  ds <- extract(ds, keep);

  ds;
});

# Keep only common samples (just in case)
sampleNames <- Reduce(intersect, lapply(dsList, FUN=getNames));
dsList <- lapply(dsList, FUN=extract, sampleNames);
print(dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Attach genotype calls to each of the above samples
#
# TO DO: There should be a general function in R.filesets that
# does this for you.  This is almost like having a 
# GenericDataFileDataFrame. /HB 2009-12-18
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (kk in seq(along=dsList)) {
  ds <- dsList[[kk]];
  tags <- setdiff(getTags(ds)[-1], c("TBN", "NGC"));
  tags <- c(tags, "NGC");
  gc <- AromaUnitGenotypeCallSet$byName(dataSet, tags=tags, chipType="*");
  setFullNamesTranslator(gc, tcgaSampleFnt);
  gc <- extract(gc, sampleNames);

  for (ii in seq(length=nbrOfFiles(ds))) {
    df <- getFile(ds, ii);
    gf <- getFile(gc, indexOf(gc, getName(df)));
    df$callFile <- gf;
  }

  rm(gc);
} # for (kk ...)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Attach total CNs to each of the above samples
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (kk in seq(along=dsList)) {
  ds <- dsList[[kk]];
  tags <- setdiff(getTags(ds)[-1], c("TBN", "NGC"));
  dsC <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType="*", paths="rawCnData/");
  setFullNamesTranslator(dsC, tcgaSampleFnt);
  dsC <- extract(dsC, sampleNames);

  for (ii in seq(length=nbrOfFiles(ds))) {
    df <- getFile(ds, ii);
    dfC <- getFile(dsC, indexOf(dsC, getName(df)));
    df$tcnFile <- dfC;
  }

  rm(dfC);
} # for (kk ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Identify all heterozygous SNPs in the genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chromosomes <- 1:24;
unitsList <- list();
for (kk in seq(along=dsList)) {
  ds <- dsList[[kk]];
  df <- getFile(ds, 1);
  gf <- df$callFile;
  
  # Identify all homozygous SNPs on this chromosome
  ugp <- getAromaUgpFile(gf);
  isHet <- isHeterozygous(gf);
  units <- whichVector(isHet);

  unitsList[[kk]] <- units;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract BAF for those SNPs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
bafList <- list();
tcnList <- list();
for (kk in seq(along=dsList)) {
  ds <- dsList[[kk]];
  df <- getFile(ds, 1);
  gf <- df$callFile;
  cf <- df$tcnFile;
  
  units <- unitsList[[kk]];

  bafs <- list();
  tcns <- list();
  for (chr in chromosomes) {
    chrTag <- sprintf("chr%02d", chr); 
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", chr, chrTag, length(chromosomes)));
    baf <- extractRawAlleleBFractions(df, chromosome=chr, units=units);    
    bafs[[chrTag]] <- baf;

    tcn <- extractRawCopyNumbers(cf, chromosome=chr);    
    tcns[[chrTag]] <- tcn;
    verbose && exit(verbose);
  }

  bafList[[kk]] <- bafs;
  tcnList[[kk]] <- tcns;
} # for (kk ...)


# Extract LH signals
lhList <- lapply(bafList, function(bafs) {
  lapply(bafs, FUN=function(baf) {
    lh <- extractRawMirroredAlleleBFractions(baf);
    # Fix new definition of LH
    lh$y <- 2*lh$y;
    lh;
  });
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Estimate smooth LHs at common loci along the genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
fitUgp <- AromaUgpFile$byChipType("GenericHuman", tags="200kb");
print(fitUgp);


lhSList <- list();
tcnSList <- list();
for (kk in seq(along=lhList)) {
  lhs <- lhList[[kk]];
  tcns <- tcnList[[kk]];

  lhSs <- list();
  tcnSs <- list();
  for (cc in seq(along=lhs)) {
    lh <- lhs[[cc]];
    tcn <- tcns[[cc]];

    chr <- lh$chromosome;
    chrTag <- sprintf("chr%02d", chr); 
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", chr, chrTag, length(lhs)));

    # Target loci
    units <- getUnitsOnChromosome(fitUgp, chromosome=chr);
    xOut <- getPositions(fitUgp, units=units);

    # Smooth LHs
    lhS <- binnedSmoothing(lh, xOut=xOut, verbose=less(verbose,5));
    lhSs[[chrTag]] <- lhS;

    # Smooth TCNs
    tcnS <- binnedSmoothing(tcn, xOut=xOut, verbose=less(verbose,5));
    tcnSs[[chrTag]] <- tcnS;

    verbose && exit(verbose);
  }

  lhSList[[kk]] <- lhSs;
  tcnSList[[kk]] <- tcnSs;
} # for (kk ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract all smoothed LH
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
L <- sapply(lhSList, FUN=function(gs) {
  s <- sapply(gs[1:23], FUN=getSignals);
  unlist(s, use.names=FALSE);
});
str(L);

T <- sapply(tcnSList, FUN=function(gs) {
  s <- sapply(gs[1:23], FUN=getSignals);
  unlist(s, use.names=FALSE);
});
T <- 2*2^T;
str(T);

CC <- sapply(lhSList, FUN=function(lhSs) {
  CC <- sapply(lhSs, FUN=function(lh) {
    rep(lh$chromosome, times=nbrOfLoci(lh));
  });
  unlist(CC, use.names=FALSE);
});
stopifnot(identical(CC[,1], CC[,2]));
CC <- CC[,1];
str(CC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Multi-source LH normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
fit <- fitPrincipalCurve(L);
print(fit);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot smoothed LH signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sampleName <- lhSList[[1]][[1]]$name;
sampleName <- gsub(",.*", "", sampleName);
lhLim <- c(0,1);


filename <- sprintf("%s,LH,raw,tracks.png", sampleName);
pathname <- file.path(figPath, filename);
if (!isFile(pathname)) {
  devNew("png", pathname, width=640, height=640);
  plotLHTracks(lhSList);
  devDone();
}


filename <- sprintf("%s,TCN,raw,tracks.png", sampleName);
pathname <- file.path(figPath, filename);
if (!isFile(pathname)) {
  devNew("png", pathname, width=640, height=640);
  plotLHTracks(tcnSList, lhLim=c(0,4), g=function(y) 2*2^y, h=2);
  devDone();
}




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# LH vs LH
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Color each chromosome differently
lhSs <- lhSList[[1]];
cols <- lapply(seq(along=lhSs), FUN=function(kk) {
  rep(kk, times=nbrOfLoci(lhSs[[kk]]));
});
cols <- unlist(cols, use.names=FALSE);

lhLabs <- unlist(tagsList);
lhLabs <- lapply(lhLabs, FUN=function(x) substitute(LH[key], list(key=as.name(x))));

filename <- sprintf("%s,LHvsLH,raw.png", sampleName);
pathname <- file.path(figPath, filename);
devNew("png", pathname, width=640, height=640);
plot(L, pch=".", col=cols, cex=2, xlim=lhLim, ylim=lhLim, xlab=lhLabs[[1]], ylab=lhLabs[[2]]);
abline(a=0, b=1, lty=2);
stext(side=3, pos=0, sampleName);
devDone();


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# TCN vs TCN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Color each chromosome differently
tcnLabs <- unlist(tagsList);
tcnLabs <- lapply(tcnLabs, FUN=function(x) substitute(TCN[key], list(key=as.name(x))));
tcnLim <- c(0,4);
	
filename <- sprintf("%s,TCNvsTCN,raw.png", sampleName);
pathname <- file.path(figPath, filename);
devNew("png", pathname, width=640, height=640);
plot(T, pch=".", col=cols, cex=2, xlim=tcnLim, ylim=tcnLim, xlab=tcnLabs[[1]], ylab=tcnLabs[[2]]);
abline(a=0, b=1, lty=2);
stext(side=3, pos=0, sampleName);
devDone();



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Normalized LHs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
LN <- backtransformPrincipalCurve(L, fit=fit, targetDimension=1);
filename <- sprintf("%s,LHvsLH,norm,robust.png", sampleName);
pathname <- file.path(figPath, filename);
devNew("png", pathname, width=640, height=640);
plot(LN, pch=".", col=cols, cex=2, xlim=lhLim, ylim=lhLim, xlab=lhLabs[[1]], ylab=lhLabs[[2]]);
abline(a=0, b=1, lty=2);
stext(side=3, pos=0, sampleName);
devDone();


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C <- T*L;

cLim <- c(0,1);
filename <- sprintf("%s,CvsC,raw.png", sampleName);
pathname <- file.path(figPath, filename);
devNew("png", pathname, width=640, height=640);
plot(C, pch=".", col=cols, cex=4, xlim=cLim, ylim=cLim, xlab="T*L", ylab="T*L");
abline(a=0, b=1, lty=2);
stext(side=3, pos=0, sampleName);
devDone();


stop()


smootherFcn <- function(lambda, xj, ..., df=5) {
   o <- order(lambda);
   lambda <- lambda[o];
   xj <- xj[o];
   fit <- robustSmoothSpline(lambda, xj, ..., df=df);
   predict(fit, x = lambda)$y;
}

fit2 <- fitPrincipalCurve(L, smoother=smootherFcn);
YN <- backtransformPrincipalCurve(L, fit=fit2, targetDimension=1);

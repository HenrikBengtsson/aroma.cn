if (interactive()) savehistory();
library("aroma.cn");
library("DNAcopy");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
fig <- 1;

figPath <- "msSeg/figures/col/";
figPath <- Arguments$getWritablePath(figPath);
figForce <- FALSE;
width <- 8*512;
height <- 400;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Chromosomes of interest
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chromosomes <- 1:25;
path <- Arguments$getReadablePath("annotationData/genomes/Human/");
pathnames <- affxparser::findFiles(path=path);
stopifnot(length(pathnames) > 0);
pathname <- pathnames[1];
db <- TabularTextFile(pathname);
print(db);
df <- readDataFrame(db);
chrLengths <- df$nbrOfBases;
names(chrLengths) <- sprintf("Chr%02d", chromosomes);
stopifnot(all(is.finite(chrLengths)));

widthScales <- chrLengths/chrLengths[1];
t <- sprintf("%.2f", widthScales);
names(t) <- names(widthScales);
print(t);

widths <- width*widthScales;
widths[widths < 420] <- 420;


if (interactive()) {
  ans <- textMenu(sprintf("Chr%02d", chromosomes));
  chromosome <- chromosomes[ans];
} else {
  chromosome <- 20;
}
chromosomes <- chromosome;

tagsList <- list("MSKCC", "Harvard", "Stanford", "Broad");
tagsList <- lapply(tagsList, FUN=c, "mscn");
dataSet <- "TCGA,GBM,testSet,pairs";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load raw CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsList <- lapply(tagsList, FUN=function(tags) {
  AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType="*");
});
dsList <- GenericDataFileSetList(dsList);
print(dsList);

sampleNames <- unlist(lapply(dsList, FUN=getNames), use.names=TRUE);
sampleNames <- sort(unique(sampleNames));
sampleNames <- "TCGA-02-0104-01Avs10A";
cat("All sample names:\n");
print(sampleNames);

tags <- Reduce(intersect, lapply(dsList, FUN=getTags));
sites <- sapply(dsList, FUN=function(ds) setdiff(getTags(ds), tags));
platforms <- sapply(dsList, FUN=getPlatform);
chipTypes <- sapply(dsList, FUN=getChipType);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment each sample and chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
methods <- list(
  HaarSeg=function(...) {
    segmentByHaarSeg(..., breaksFdrQ=1e-6)
  }, 
  CBS=function(...) {
    segmentByCBS(..., weights=NULL)
  },
  wCBS=function(...) {
    segmentByCBS(...)
  }
);

segCols <- c(HaarSeg=2, CBS=3, wCBS=4, GLAD=5);

srcCols <- gray(seq(from=0.3, to=0.8, length.out=5));
xScale <- 1; # xScale <- 1e-6;	

  


for (ii in seq(along=sampleNames)) {
  sampleName <- sampleNames[ii];
  verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                                      ii, sampleName, length(sampleNames)));
  
  dfList <- getFileList(dsList, name=sampleName);
  print(dfList);
  
  for (cc in seq(along=chromosomes)) {
    chrTag <- sprintf("Chr%02d", chromosome);
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", 
                                         cc, chrTag, length(chromosomes)));
    chromosome <- chromosomes[cc];
  
    cnM <- extractMergedRawCopyNumbers(dfList, chromosome=chromosome, 
                                             bandwidth=200e3, verbose=verbose);
    verbose & print(verbose, cnM);
    nbrOfSources <- length(unique(cnM$id));
    
    figName <- sprintf("%s,%s,segmentation", sampleName, chrTag);
    filename <- sprintf("%s.png", figName);
    pathname <- Arguments$getWritablePathname(filename, path=figPath, 
                                                      mustNotExist=FALSE);
    if (figForce || !isFile(pathname)) {
#      devSet(figName, width=8, height=4);
      devNew(png, pathname, width=widths[chromosome], height=height);

      par(mar=c(3.5,3.5,2,1)+0.1, mgp=c(2.3,0.7,0), cex.lab=2, cex.axis=2);
    
      plot(cnM, col=srcCols[cnM$id], pch=".", cex=2, xScale=xScale);
      cex <- 0.8;
      cex <- 2;
      stext(side=3, pos=0, cex=2, sampleName);
      label <- sprintf("Chr%02d (n=%d)", chromosome, nbrOfLoci(cnM));
      stext(side=3, pos=1, cex=2, label);
      stext(side=3, pos=1, line=-1, cex=cex, sprintf("s*=%.3f", cnM$sigma));
      legend("topleft", lwd=3, col=segCols, names(methods), horiz=TRUE, 
             bty="n", x.intersp=0.5, cex=cex);
      
      fitList <- list();
      for (mm in seq(along=methods)) {
        methodName <- names(methods)[mm];
        segmentByFcn <- methods[[mm]];
        verbose && enter(verbose, sprintf("Segmentation #%d ('%s') of %d", 
                                             mm, methodName, length(methods)));
        fit <- segmentByFcn(cnM, verbose=less(verbose, 5));
        fit$call <- NULL;
        fitList[[methodName]] <- fit;
        drawCnRegions(fit, lwd=5, col="white", xScale=xScale);
        drawCnRegions(fit, lwd=3, col=segCols[methodName], xScale=xScale);
        t <- attr(fit, "processingTime");
        verbose && cat(verbose, "Processing time:");
        verbose && print(verbose, t);
        verbose && cat(verbose, "Processing time per 1,000 loci:");
        verbose && print(verbose, 1000*t/nbrOfLoci(cnM));
        verbose && exit(verbose);
      } # for (mm ...)
      
      legend("topleft", lwd=3, col=segCols, names(methods), horiz=TRUE, 
             bty="n", x.intersp=0.5, cex=cex);
      
      devDone();
    } # if (!figForce && isFile(pathname))

    verbose && exit(verbose);
    } # for (cc ...)

  verbose && exit(verbose);
} # for (ii ...)

############################################################################
# HISTORY:
# 2009-05-13
# o Added generation of PNGs.
# 2009-05-11
# o Created.
############################################################################ 

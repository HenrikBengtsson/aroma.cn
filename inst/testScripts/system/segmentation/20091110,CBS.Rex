if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");

figPath <- Arguments$getWritablePath("figures");

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
# Load CN data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*");
print(ds);

ugp <- getAromaUgpFile(ds);
print(ugp);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Segment the total CNs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Clim <- c(0,5);
Clab <- expression(C);
xScale <- 1e-6;

chromosomes <- getChromosomes(ugp);

for (ii in seq(ds)) {
  df <- getFile(ds, ii);
  sampleName <- getName(df);
  verbose && enter(verbose, sprintf("Sample %d ('%s') of %d", 
                            ii, sampleName, length(ds)));

  # Extract all CN data
  C <- df[,1,drop=TRUE];

  # Transform?
  if (!hasTag(df, "logRatio")) {
    C <- 2*10^C;
    Clim <- c(0,5);
  } else {
    Clim <- c(-3,3);
  }

  x <- getPositions(ugp);

  for (cc in seq(along=chromosomes)) {
    chromosome <- chromosomes[cc];
    chrTag <- sprintf("Chr%02d", chromosome);
    verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", 
                              cc, chrTag, length(chromosomes)));
  
    # Identify units
    units <- getUnitsOnChromosome(ugp, chromosome);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    # To few data points?
    if (length(units) < 200) {
      verbose && cat(verbose, "Too few data points: ", length(units));
      verbose && exit(verbose);
      next;
    }
  
    # Get genomic positions
    cn <- RawCopyNumbers(x=x[units], cn=C[units], chromosome=chromosome, name=sampleName);
    plot(cn, ylim=Clim, xScale=xScale);
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, chrTag);

    # Segmentation
    verbose && enter(verbose, "Segmenting");
    fit <- segmentByCBS(cn, cache=TRUE);
    verbose && print(verbose, fit);
    cnr <- extractCopyNumberRegions(fit);
    verbose && print(verbose, cnr);
    verbose && exit(verbose);

    # Plot
    drawLevels(cnr, lwd=4, xScale=xScale);

    verbose && exit(verbose);
  } # for (cc ...)

  verbose && exit(verbose);
} # for (ii ...)

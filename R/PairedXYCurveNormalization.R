###########################################################################/**
# @RdocClass PairedXYCurveNormalization
#
# @title "The PairedXYCurveNormalization class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#  \item{ds}{An @see "aroma.core::AromaUnitTotalCnBinarySet" of 
#     "test" samples to be normalized.}
#  \item{dsTarget}{An @see "aroma.core::AromaUnitTotalCnBinarySet" of 
#     paired target samples.}
#  \item{subsetToFit}{The subset of loci to be used to fit the 
#    normalization functions.
#    If @NULL, loci on chromosomes 1-22 are used, but not on ChrX and ChrY.
#  }
#  \item{flavor}{A @character string specifying the type of 
#     correction applied.}
#  \item{tags}{(Optional) Sets the tags for the output data sets.}
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{Henrik Bengtsson}
#*/########################################################################### 
setConstructorS3("PairedXYCurveNormalization", function(ds=NULL, dsTarget=NULL, subsetToFit=NULL, flavor=c("v1"), tags="*", ...) {
  # Validate arguments
  if (!is.null(ds)) {
    # Argument 'flavor':
    flavor <- match.arg(flavor);

    # Arguments 'ds' and 'dsTarget'
    dsList <- list(ds=ds, dsTarget=dsTarget);
    className <- "AromaUnitTotalCnBinarySet";
    for (kk in seq(along=dsList)) {
      key <- names(dsList)[kk];
      ds <- dsList[[kk]];
      if (!inherits(ds, className)) {
        throw(sprintf("Argument '%s' is not of class %s: %s", key, 
                                         className, class(ds)[1]));
      }
    } # for (kk ...)

    # Assert that each data set contains the same number of files
    for (jj in 1:(length(dsList)-1)) {
      keyJJ <- names(dsList)[jj];
      dsJJ <- dsList[[jj]];
      nJJ <- nbrOfFiles(dsJJ);
      chipTypeJJ <- getChipType(dsJJ);
      for (kk in (jj+1):length(dsList)) {
        keyKK <- names(dsList)[kk];
        dsKK <- dsList[[kk]];
        nKK <- nbrOfFiles(dsKK);
        chipTypeKK <- getChipType(dsKK);

        # Assert that each data set contains the same number of files
        if (nKK != nJJ) {
          throw(sprintf("The number of files in '%s' and '%s' does not match: %s != %s", keyKK, keyJJ, nKK, nJJ));
        }

        # Assert that each data set is for the same chip type
        if (chipTypeKK != chipTypeJJ) {
          throw(sprintf("The chip types for '%s' and '%s' does not match: %s != %s", keyKK, keyJJ, chipTypeKK, chipTypeJJ));
        }
      } # for (kk ...)
    } # for (jj ...)

    # Assert that the UGP file exists
    ugp <- getAromaUgpFile(ds);

    # Argument 'subsetToFit':
    if (is.null(subsetToFit)) {
    } else if (is.character(subsetToFit)) {
      throw("Yet not implemented: Argument 'subsetToFit' is of type character.");
    } else {
      subsetToFit <- Arguments$getIndices(subsetToFit, 
                                          range=c(1, nbrOfUnits(ugp)));
    }
  } # if (!is.null(ds))

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  this <- extend(Object(...), "PairedXYCurveNormalization",
    .ds = ds,
    .dsTarget = dsTarget,
    .subsetToFit = subsetToFit,
    .flavor = flavor
  );

  setTags(this, tags);

  this;
})


setMethodS3("as.character", "PairedXYCurveNormalization", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  s <- c(s, sprintf("Flavor: %s", getFlavor(this)));

  dsList <- getDataSets(this);

  dsList <- getDataSets(this);
  s <- c(s, sprintf("Data sets (%d):", length(dsList)));
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    s <- c(s, sprintf("<%s>:", capitalize(names(dsList)[kk])));
    s <- c(s, as.character(ds));
  } 
 
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getAsteriskTags", "PairedXYCurveNormalization", function(this, collapse=NULL, ...) {
  tags <- "PXYCN";

  flavor <- getFlavor(this);
  if (flavor != "v1") {
    tags <- c(tags, flavor);
  }

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }
  
  tags;
}, private=TRUE)


setMethodS3("getName", "PairedXYCurveNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getSubsetToFit", "PairedXYCurveNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  units <- this$.subsetToFit;
  if (is.null(units)) {
    verbose && enter(verbose, "Identify subset of units for fitting the normalization function");

    verbose && enter(verbose, "Retrieving the UGP file");
    ds <- getInputDataSet(this);
    ugp <- getAromaUgpFile(ds);
    verbose && print(verbose, ugp);
    verbose && exit(verbose);
  
    verbose && enter(verbose, "Querying UGP for units on chromosomes of interest");
    chromosomes <- 1:22;
    verbose && cat(verbose, "Chromosomes to fit: ", 
                                             seqToHumanReadable(chromosomes));
    units <- sapply(chromosomes, FUN=function(cc) {
      getUnitsOnChromosome(ugp, cc);
    });
    units <- unlist(units, use.names=FALSE);
    units <- unique(units);
    units <- sort(units);
    verbose && str(verbose, units);
    verbose && exit(verbose);

    this$.subsetToFit <- units;

    verbose && exit(verbose);
  }

  units;
}, protected=TRUE)


setMethodS3("getFlavor", "PairedXYCurveNormalization", function(this, ...) {
  this$.flavor;
}, protected=TRUE)


setMethodS3("getTags", "PairedXYCurveNormalization", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  ds <- getInputDataSet(this);
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("setTags", "PairedXYCurveNormalization", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})

 
setMethodS3("getFullName", "PairedXYCurveNormalization", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getDataSets", "PairedXYCurveNormalization", function(this, ...) {
  list(test=this$.ds, target=this$.dsTarget);
}, protected=TRUE)

setMethodS3("getInputDataSet", "PairedXYCurveNormalization", function(this, ...) {
  this$.ds;
})

setMethodS3("getTargetDataSet", "PairedXYCurveNormalization", function(this, ...) {
  this$.dsTarget;
})

setMethodS3("getRootPath", "PairedXYCurveNormalization", function(this, ...) {
  "totalAndFracBData";
})

setMethodS3("getPath", "PairedXYCurveNormalization", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ds <- getInputDataSet(this);
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");

  # Verify that it is not the same as the input path
  inPath <- getPath(getInputDataSet(this));
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  # Create path?
  if (create) {
    if (!isDirectory(path)) {
      mkdirs(path);
      if (!isDirectory(path))
        throw("Failed to create output directory: ", path);
    }
  }

  path;
})


setMethodS3("nbrOfFiles", "PairedXYCurveNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  nbrOfFiles(ds);
})


setMethodS3("getOutputDataSet", "PairedXYCurveNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- fromFiles(ds, path=path, ...);
  res;
})


setMethodS3("getPairedDataSet", "PairedXYCurveNormalization", function(this, array, ..., verbose=FALSE) {
  ds <- getInputDataSet(this);
  nbrOfArrays <-nbrOfArrays(ds);

  # Argument 'array':
  array <- Arguments$getIndex(array, range=c(1, nbrOfArrays));

  df <- getFile(ds, array);
  name <- getName(df);

  verbose && enter(verbose, sprintf("Extracting paired data set for array %d ('%s') of %d", array, name, nbrOfArrays));
  dsT <- getTargetDataSet(this);
  dfT <- getFile(dsT, array);

  dsPair <- newInstance(this, list(df, dfT));
  verbose && cat(verbose, "Pair:");
  verbose && print(verbose, dsPair);

  verbose && exit(verbose);

  dsPair;
}, protected=TRUE)



setMethodS3("process", "PairedXYCurveNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Paired (x,y)-curve normalization");
  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of arrays: ", nbrOfFiles);

  flavor <- getFlavor(this);
  verbose && cat(verbose, "Flavor: ", flavor);

  ds <- getInputDataSet(this);
  chipType <- getChipType(ds, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);

  outPath <- getPath(this);

  units <- NULL;
  for (kk in seq(length=nbrOfFiles)) {
    dsPair <- getPairedDataSet(this, array=kk, verbose=less(verbose,5));
    verbose && cat(verbose, "Pair:");
    verbose && print(verbose, dsPair);

    df <- getFile(dsPair, 1);
    name <- getName(df);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                                                    kk, name, nbrOfFiles));

    # Output file
    filename <- getFilename(df);
    pathname <- Arguments$getReadablePathname(filename, path=outPath, mustExist=FALSE);

    # Nothing to do?
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already normalized.");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Reading all data");
    theta <- extractMatrix(dsPair, verbose=less(verbose,5));
    verbose && str(verbose, theta);
    verbose && exit(verbose);
    nbrOfUnits <- nrow(theta);

    subsetToFit <- getSubsetToFit(this);
    verbose && cat(verbose, "Subset to fit:");
    verbose && str(verbose, subsetToFit);

    verbose && enter(verbose, "Fitting");
    thetaFit <- theta[subsetToFit,,drop=FALSE];
    fit <- fitPrincipalCurve(thetaFit, ..., verbose=less(verbose,10));
    rm(thetaFit, subsetToFit);
    verbose && str(verbose, fit);
    verbose && exit(verbose);

    verbose && enter(verbose, "Backtransforming data (normalizing)");
    thetaN <- backtransformPrincipalCurve(theta, fit=fit, 
                         targetDimension=2, verbose=less(verbose, 1));
    thetaN <- thetaN[,1,drop=TRUE];
    rm(fit);
    verbose && str(verbose, thetaN);
    verbose && exit(verbose);

    # Sanity check
    stopifnot(length(thetaN) == nbrOfUnits);

    verbose && enter(verbose, "Storing normalized data");

    verbose && enter(verbose, "Allocating to temporary file");
    pathnameT <- sprintf("%s.tmp", pathname);
    outPath <- Arguments$getWritablePath(outPath);
    file.copy(getPathname(df), pathnameT);
    dfN <- newInstance(df, pathnameT);
    srcFiles <- lapply(dsPair, function(df) {
      list(
        filename = getFilename(df),
        filesize = getFileSize(df),
        checksum = getChecksum(df)
      )
    });
    footer <- readFooter(dfN);
    footer$srcFiles <- srcFiles;
    writeFooter(dfN, footer);
    rm(srcFiles, footer);
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing to temporary file");
    dfN[units,1] <- thetaN;
    rm(thetaN);
    verbose && exit(verbose);

    # Renaming
    verbose && enter(verbose, "Renaming temporary file");
    file.rename(pathnameT, pathname);
    if (!isFile(pathname)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    verbose && exit(verbose);
    verbose && exit(verbose);

    # More clean up
    rm(df, dfN, pathnameT);

    verbose && exit(verbose);
  } # for (kk ...)

  res <- getOutputDataSet(this, verbose=less(verbose, 1)); 

  verbose && exit(verbose);

  invisible(res);
})

############################################################################
# HISTORY:
# 2009-07-15
# o Created.
############################################################################ 

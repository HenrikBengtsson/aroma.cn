###########################################################################/**
# @RdocClass TumorBoostNormalization
#
# @title "The TumorBoostNormalization class"
#
# \description{
#  @classhierarchy
#
#  TumorBoost is normalization method that normalizes the allele B fractions
#  of a tumor sample given the allele B fractions and genotype calls for
#  a matched normal.
#  The method is a single-sample (single-pair) method.  It does not require
#  total copy number estimates.  
#  The normalization is done such that the total copy number is unchanged
#  afterwards.
# }
# 
# @synopsis
#
# \arguments{
#  \item{dsT}{...}
#  \item{dsN}{...}
#  \item{gcN}{...}
#  \item{flavor}{Type of correction applied..}
#  \item{tags}{(Optional) Sets the tags for the output data sets.}
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("TumorBoostNormalization", function(dsT=NULL, dsN=NULL, gcN=NULL, flavor=c("v1", "v2"), tags="*", ...) {
  # Validate arguments
  if (!is.null(dsT)) {
    # Argument 'flavor':
    flavor <- match.arg(flavor);

    # Arguments 'dsT' and 'dsN'
    dsList <- list(dsT=dsT, dsN=dsN);
    className <- "AromaUnitFracBCnBinarySet";
    for (kk in seq(along=dsList)) {
      key <- names(dsList)[kk];
      ds <- dsList[[kk]];
      if (!inherits(ds, className)) {
        throw(sprintf("Argument '%s' is not of class %s: %s", key, 
                                         className, class(ds)[1]));
      }
    }

    # Argument 'gcN':
    className <- "AromaUnitGenotypeCallSet";
    if (!inherits(gcN, className)) {
      throw(sprintf("Argument '%s' is not of class %s: %s", key, 
                                         className, class(gcN)[1]));
    }

    # Assert that each data set contains the same number of files
    dsList$gcN <- gcN;
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
  } # if (!is.null(dsT))

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  this <- extend(Object(...), "TumorBoostNormalization",
    .dsT = dsT,
    .dsN = dsN,
    .gcN = gcN,
    .flavor = flavor
  );

  setTags(this, tags);

  this;
})


setMethodS3("as.character", "TumorBoostNormalization", function(x, ...) {
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


setMethodS3("getAsteriskTags", "TumorBoostNormalization", function(this, collapse=NULL, ...) {
  tags <- "TBN";

  flavor <- getFlavor(this);
  if (flavor != "v1") {
    tags <- c(tags, flavor);
  }

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }
  
  tags;
}, private=TRUE)


setMethodS3("getName", "TumorBoostNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getFlavor", "TumorBoostNormalization", function(this, ...) {
  this$.flavor;
}, protected=TRUE)


setMethodS3("getTags", "TumorBoostNormalization", function(this, collapse=NULL, ...) {
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


setMethodS3("setTags", "TumorBoostNormalization", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})

 
setMethodS3("getFullName", "TumorBoostNormalization", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getDataSets", "TumorBoostNormalization", function(this, ...) {
  list(tumor=this$.dsT, normal=this$.dsN, normalCalls=this$.gcN);
}, protected=TRUE)

setMethodS3("getInputDataSet", "TumorBoostNormalization", function(this, ...) {
  this$.dsT;
})

setMethodS3("getNormalDataSet", "TumorBoostNormalization", function(this, ...) {
  this$.dsN;
})

setMethodS3("getNormalGenotypeCallSet", "TumorBoostNormalization", function(this, ...) {
  this$.gcN;
})

setMethodS3("getRootPath", "TumorBoostNormalization", function(this, ...) {
  "totalAndFracBData";
})

setMethodS3("getPath", "TumorBoostNormalization", function(this, create=TRUE, ...) {
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


setMethodS3("nbrOfFiles", "TumorBoostNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  nbrOfFiles(ds);
})


setMethodS3("getOutputDataSet", "TumorBoostNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- fromFiles(ds, path=path, ...);
  res;
})


setMethodS3("process", "TumorBoostNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  units <- NULL;

  verbose && enter(verbose, "TumorBoost normalization");
  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of arrays: ", nbrOfFiles);

  flavor <- getFlavor(this);
  verbose && cat(verbose, "Flavor: ", flavor);

  outPath <- getPath(this);
  dsList <- getDataSets(this);
  for (kk in seq(length=nbrOfFiles)) {
    dfList <- lapply(dsList, FUN=getFile, kk);
    dfT <- dfList$tumor;
    name <- getName(dfT);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                                                    kk, name, nbrOfFiles));

    # Output file
    filename <- getFilename(dfT);
    pathname <- Arguments$getReadablePathname(filename, path=outPath, mustExist=FALSE);

    # Nothing to do?
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already normalized.");
      verbose && exit(verbose);
      next;
    }

    verbose && print(verbose, dfList);

    if (is.null(units)) {
      verbose && enter(verbose, "Identifying units to read");
      units <- seq(length=nbrOfUnits(dfList$tumor));
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Reading all data");
    betaT <- dfList$tumor[units,1,drop=TRUE];
    keep <- is.finite(betaT);
    unitsT <- units[keep];
    betaT <- betaT[keep];
    verbose && cat(verbose, "Allele B fractions for the tumor:");
    verbose && str(verbose, betaT);

    verbose && cat(verbose, "Allele B fractions for the matched normal:");
    betaN <- dfList$normal[unitsT,1,drop=TRUE];
    verbose && str(verbose, betaN);

    verbose && cat(verbose, "Genotypes for the matched normal:");
    gfN <- dfList$normalCalls;
    muN <- extractGenotypes(gfN, units=unitsT, encoding="fracB", drop=TRUE);
    verbose && str(verbose, muN);
    verbose && exit(verbose);


    verbose && enter(verbose, "Normalizing tumor allele B fractions");
    verbose && cat(verbose, "Flavor: ", flavor);
    verbose && enter(verbose, "Estimating SNP effects");
    delta <- (betaN - muN);
    verbose && str(verbose, delta);
    verbose && exit(verbose);

    verbose && enter(verbose, "Rescaling correction factor");
    if (flavor == "v1") {
      b <- 1;
    } else if (flavor == "v2") {
      b <- rep(1, length(delta));
      isUp <- (betaT > betaN);
      idxs <- whichVector(isUp);
      b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
      idxs <- whichVector(!isUp);
      b[idxs] <- betaT[idxs]/betaN[idxs];
      rm(isUp,idxs);
    }
    verbose && cat(verbose, "Scaling factor:");
    verbose && str(verbose, b);
    verbose && summary(verbose, b);
    verbose && exit(verbose);

    verbose && enter(verbose, "Normalizing");
    betaTC <- betaT - b*delta;
    verbose && str(verbose, betaTC);
    verbose && exit(verbose);
    verbose && exit(verbose);


    verbose && enter(verbose, "Storing normalized data");

    verbose && enter(verbose, "Allocating to temporary file");
    pathnameT <- sprintf("%s.tmp", pathname);
    dfT <- dfList$tumor;
    outPath <- Arguments$getWritablePath(outPath);
    file.copy(getPathname(dfT), pathnameT);
    dfTC <- newInstance(dfT, pathnameT);
    srcFiles <- lapply(dfList, FUN=function(df) {
      list(
        filename = getFilename(dfT),
        filesize = getFileSize(dfT),
        checksum = getChecksum(dfT)
      )
    });
    footer <- readFooter(dfTC);
    footer$srcFiles <- footer;
    writeFooter(dfTC, footer);
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing to temporary file");
    dfTC[unitsT,1] <- betaTC;
    verbose && exit(verbose);

    # Renaming
    verbose && enter(verbose, "Renaming temporary file");
    file.rename(pathnameT, pathname);
    if (!isFile(pathname)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    verbose && exit(verbose);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)

  res <- getOutputDataSet(this, verbose=less(verbose, 1)); 

  verbose && exit(verbose);

  invisible(res);
})

############################################################################
# HISTORY:
# 2009-06-22
# o Added model 'flavor' "v2".
# 2009-06-08
# o The constructor of TumorBoostNormalization now only takes an
#   AromaUnitGenotypeCallSet for argument 'gcN'.  It no longer takes an
#   AromaUnitFracBCnBinarySet object.
# 2009-05-17
# o Now the constructor of TumorBoostNormalization asserts that there are
#   no stray arguments.
# 2009-04-29
# o Created.
############################################################################ 

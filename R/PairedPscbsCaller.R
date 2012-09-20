setConstructorS3("PairedPscbsCaller", function(dataSet=NULL, calls=c("ROH", "AB", "LOH"), ...) {
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "PairedPSCBSFileSet");
  }

  # Argument 'calls':
  calls <- match.arg(calls, several.ok=TRUE);

  extend(AromaTransform(dataSet=dataSet, ...,
               .reqSetClass="PairedPSCBSFileSet"), "PairedPscbsCaller",
    .calls = calls
  );
}) # PairedPscbsCaller()



setMethodS3("getAsteriskTags", "PairedPscbsCaller", function(this, collapse=NULL, ...) {
  calls <- this$.calls;
  tags <- c("call", calls);

  # Collapsed or split?
  tags <- Arguments$getTags(tags, collapse=collapse);

  tags;
}, private=TRUE)


setMethodS3("getRootPath", "PairedPscbsCaller", function(this, ...) {
  "pscbsData";
}, protected=TRUE)

setMethodS3("getPath", "PairedPscbsCaller", function(this, create=TRUE, ...) {
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


setMethodS3("getParameters", "PairedPscbsCaller", function(this, ...) {
  calls <- this$.calls;
  params <- list(calls=calls);
  params;
}, private=TRUE)



setMethodS3("process", "PairedPscbsCaller", function(this, ..., force=FALSE) {
  # Argument 'force':
  force <- Arguments$getLogical(force);

  sms <- getInputDataSet(this);

  verbose && enter(verbose, "Calling LOH and AB");

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already done. Skipping");
    res <- getOutputDataSet(this);
    verbose && exit(verbose);
    return(res);
  }

  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, sms);

  pathD <- getPath(this);
  verbose && cat(verbose, "Output path: ", pathD);

  verbose && cat(verbose, "Number of samples: ", length(sms));
  
  for (ii in seq(sms)) {
    smf <- getFile(sms, ii);
    sampleName <- getName(smf);
    verbose && enter(verbose, "Tumor-normal pair #%d ('%s') of %d", ii, sampleName, length(sms));
  
    filename <- getFilename(smf);
    pathname <- file.path(pathD, filename);

    # Sanity check
    stopifnot(getAbsolutePath(pathname) != getAbsolutePath(getFullName(smf)));

    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already called. Skipping.");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Loading segmentation data");
    fit <- loadObject(getPathname(smf));

    # Sanity check
    fit <- Arguments$getInstanceOf(fit, "PairedPSCBS");
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Calling ROH");
    fit <- callROH(fit, verbose=less(verbose, 5));   
    verbose && exit(verbose);
  
    verbose && enter(verbose, "Calling AB");
    fit <- callAB(fit, verbose=less(verbose, 5));   
    verbose && exit(verbose);
  
    verbose && enter(verbose, "Saving");
    saveObject(fit, file=pathname);
    verbose && exit(verbose);
  
    verbose && exit(verbose);
  } # for (ii ...)

  res <- getOutputDataSet(this);
  verbose && print(verbose, res);
  verbose && exit(verbose);

  res;
}) # process()


# AD HOC
setMethodS3("getPlatform", "PairedPscbsCaller", function(this, ...) {
  "GenericPlatform";
}, protected=TRUE)



##########################################################################
# HISTORY:
# 2012-09-19
# o Made as an AromaTransform for now.
# o Created.
##########################################################################

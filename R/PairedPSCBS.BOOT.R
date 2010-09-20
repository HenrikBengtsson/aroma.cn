setMethodS3("bootstrap", "PairedPSCBS", function(fit, B=100, flavor=c("c"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(1,Inf));

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bootstrap
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (flavor == "c") {
    resampleFcn <- resampleC;
  }

  fitBList <- list();
  for (bb in 1:100) {
    print(bb);
    fitBList[[bb]] <- resampleFcn(fit);
  }

  fitBList;  
}) # bootstrap()



setMethodS3("resampleC", "PairedPSCBS", function(fit, by=c("betaTN", "betaT"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'by':
  by <- match.arg(by);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Resample (TCN,BAF) signals and reestimate segmentation means");


  data <- fit$data;
  snpFields <- c("betaN", "betaT", "betaTN", "muN");
  npFields <- "CT";
  snpFields <- intersect(names(data), snpFields);
  npFields <- intersect(names(data), npFields);
  verbose && cat(verbose, "SNP fields to be resampled:");
  verbose && print(verbose, snpFields);
  verbose && cat(verbose, "Non-polymorphic fields to be resampled:");
  verbose && print(verbose, npFields);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample TCN segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  x <- data$x;
  segs <- fit$output;
  nbrOfSegments <- nrow(segs);

  verbose && enter(verbose, "Resampling TCN segments");
  ids <- unique(segs[["tcn.id"]]);
  for (ii in seq(along=ids)) {
    id <- ids[ii];
    verbose && enter(verbose, sprintf("TCN segment #%d of %d", ii, length(ids)));

    # Identify all subsegments for this TCN segment
    idxs <- which(segs[["tcn.id"]] == id);
    n <- length(idxs);

    segsII <- segs[idxs,,drop=FALSE];

    # Identify loci in segment
    start <- segsII$tcn.loc.start[1];
    stop <- segsII$tcn.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    if (length(units) > 1) {
      # Resample
      unitsS <- sample(units, replace=TRUE);
  
      # Resample data
      for (field in npFields) {
        data[[field]][units] <- data[[field]][unitsS];
      }
    }

    verbose && exit(verbose);
  } # for (ii ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample DH segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (jj in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, nbrOfSegments));
    segsJJ <- segs[jj,,drop=FALSE];

    # Identify loci in segment
    start <- segsJJ$dh.loc.start[1];
    stop <- segsJJ$dh.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    # Identify SNPs and non-SNPs
    muN <- data$muN[units];
    isSnp <- is.finite(muN);

    # Identify heterozygous SNPs
    isHet <- (muN == 1/2);
    hets <- whichVector(isSnp & isHet);

    units <- units[hets];

    if (length(units) > 1) {
      # Resample
      unitsS <- sample(units, replace=TRUE);
  
      # Resample data
      for (field in snpFields) {
        data[[field]][units] <- data[[field]][unitsS];
      }
    }

    verbose && exit(verbose);
  } # for (jj ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the segmentation parameter estimates, e.g. mean levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  x <- data$x;

  verbose && enter(verbose, "Updating TCN fields");
  ids <- unique(segs[["tcn.id"]]);
  for (kk in seq(along=ids)) {
    id <- ids[kk];
    verbose && enter(verbose, sprintf("TCN segment #%d of %d", kk, length(ids)));

    # Identify all subsegments for this TCN segment
    idxs <- which(segs[["tcn.id"]] == id);
    n <- length(idxs);

    segsKK <- segs[idxs,,drop=FALSE];

    # Identify loci in segment
    start <- segsKK$tcn.loc.start[1];
    stop <- segsKK$tcn.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    # Update TCN mean level
    y <- data$CT[units];
    yMean <- mean(y, na.rm=TRUE);
    segs[idxs,"tcn.mean"] <- rep(yMean, times=length(idxs));

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  verbose && enter(verbose, "Updating DH fields");
  for (kk in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", kk, nbrOfSegments));

    segsKK <- segs[kk,,drop=FALSE];

    # Identify loci in segment
    start <- segsKK$dh.loc.start[1];
    stop <- segsKK$dh.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    # Identify SNPs and non-SNPs
    muN <- data$muN[units];
    isSnp <- is.finite(muN);

    # Identify heterozygous SNPs
    isHet <- (muN == 1/2);
    hets <- whichVector(isSnp & isHet);

    # Update DH mean level
    beta <- data[[by]][units];
    y <- 2*abs(beta[isHet] - 1/2);
    yMean <- mean(y, na.rm=TRUE);
    segs[kk,"dh.mean"] <- yMean;

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Store results
  fitB <- fit;
  fitB$data <- data;
  fitB$output <- segs;

  verbose && exit(verbose);

  fitB;
}) # resampleC()




setMethodS3("resampleA", "PairedPSCBS", function(fit, by=c("betaTN", "betaT"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'by':
  by <- match.arg(by);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Resample (TCN,BAF) signals and reestimate segmentation means");



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify nested TCN units and DH units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- fit$data;
  x <- data$x;
  segs <- fit$output;
  nbrOfSegments <- nrow(segs);

  verbose && enter(verbose, "Identifying nested TCN units and DH units");
  unitsList <- list();
  ids <- unique(segs[["tcn.id"]]);
  for (ii in seq(along=ids)) {
    id <- ids[ii];
    verbose && enter(verbose, sprintf("TCN segment #%d of %d", ii, length(ids)));
    unitsII <- list();

    # Identify all subsegments for this TCN segment
    idxs <- which(segs[["tcn.id"]] == id);
    n <- length(idxs);

    segsII <- segs[idxs,,drop=FALSE];

    # Identify loci in segment
    start <- segsII$tcn.loc.start[1];
    stop <- segsII$tcn.loc.end[1];
    units <- whichVector(start <= x & x <= stop);
    tcnUnits <- units;

    # For each DH segment in this TCN segment        
    dhList <- list();
    for (jj in seq(length=n)) {
      verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, n));
      segsJJ <- segsII[jj,,drop=FALSE];

      # Identify loci in segment
      start <- segsJJ$dh.loc.start[1];
      stop <- segsJJ$dh.loc.end[1];
      units <- whichVector(start <= x & x <= stop);

      # Identify SNPs and non-SNPs
      muN <- data$muN[units];
      isSnp <- is.finite(muN);
      nonSnps <- whichVector(!isSnp);

      # Identify heterozygous SNPs
      isHet <- (muN == 1/2);
      hets <- whichVector(isSnp &  isHet);
      homs <- whichVector(isSnp & !isHet);

      dhList[[jj]] <- list(dhUnits=units, nonSnps=units[nonSnps], hets=units[hets], homs=units[homs]);
      verbose && exit(verbose);
    } # for (jj ...)
    dhHetsList <- lapply(dhList, FUN=function(x) x$hets);

    unitsList[[ii]] <- list(tcnUnits=tcnUnits, dhHetsList=dhHetsList);

    verbose && exit(verbose);
  } # for (ii ...)

  verbose && cat(verbose, "Units per segment:");
  verbose && str(verbose, unitsList);
  verbose && cat(verbose, "All units in TCN segment:");
  units <- unlist(lapply(unitsList, FUN=function(x) x[[1]]), use.names=FALSE);
  units <- sort(units);
  verbose && str(verbose, units);
  verbose && cat(verbose, "All heterozygous SNPs in DH segment:");
  units <- unlist(lapply(unitsList, FUN=function(x) x[[2]]), use.names=FALSE);
  units <- sort(units);
  verbose && str(verbose, units);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample locus indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Resampling unit indices");
  unitsListS <- list();
  for (ii in seq(along=unitsList)) {
    verbose && enter(verbose, sprintf("TCN segment #%d of %d", ii, length(unitsList)));

    unitsListII <- unitsList[[ii]];
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, unitsListII);

    # (a) TCN units
    tcnUnits <- unitsListII$tcnUnits;

    # (b) Heterozygous SNPs
    dhHetsList <- unitsListII$dhHetsList;

    # (c) TCN units that are non heterozygous SNPs
    dhHets <- unlist(dhHetsList, use.names=FALSE);
    nonDhHets <- setdiff(tcnUnits, dhHets);

    # Resample (b)
    dhHetsList <- lapply(dhHetsList, FUN=sample, replace=TRUE);

    # Resample (c)
    nonDhHets <- sample(nonDhHets, replace=TRUE);

    # "Resample" (a)
    tcnUnits <- c(nonDhHets, dhHets);

    # Stop
    unitsListII$dhHetsList <- dhHetsList;
    unitsListII$tcnUnits <- tcnUnits;

    unitsListS[[ii]] <- unitsListII;
    verbose && exit(verbose);
  } # for (ii ...)

  verbose && cat(verbose, "Units per segment:");
  verbose && str(verbose, unitsListS);
  verbose && cat(verbose, "All units in TCN segment:");
  units <- unlist(lapply(unitsListS, FUN=function(x) x[[1]]), use.names=FALSE);
  units <- sort(units);
  verbose && str(verbose, units);
  verbose && cat(verbose, "All heterozygous SNPs in DH segment:");
  units <- unlist(lapply(unitsListS, FUN=function(x) x[[2]]), use.names=FALSE);
  units <- sort(units);
  verbose && str(verbose, units);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Resampling data");

  units <- unlist(lapply(unitsList, FUN=function(x) x[[1]]), use.names=FALSE);
  unitsS <- unlist(lapply(unitsListS, FUN=function(x) x[[1]]), use.names=FALSE);
  idxs <- seq(along=units);
  idxs[units] <- unitsS;

  fields <- c("CT", "betaN", "betaT", "betaTN", "muN");
  fields <- intersect(names(data), fields);
  verbose && cat(verbose, "Fields:");
  verbose && print(verbose, fields);
  for (ff in seq(along=fields)) {
    data[[ff]] <- data[[ff]][idxs];
  }
  data$units <- idxs;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the segmentation parameter estimates, e.g. mean levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  x <- data$x;

  verbose && enter(verbose, "Updating TCN fields");
  ids <- unique(segs[["tcn.id"]]);
  for (kk in seq(along=ids)) {
    id <- ids[kk];
    verbose && enter(verbose, sprintf("TCN segment #%d of %d", kk, length(ids)));

    # Identify all subsegments for this TCN segment
    idxs <- which(segs[["tcn.id"]] == id);
    n <- length(idxs);

    segsKK <- segs[idxs,,drop=FALSE];

    # Identify loci in segment
    start <- segsKK$tcn.loc.start[1];
    stop <- segsKK$tcn.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    # Update TCN mean level
    y <- data$CT[units];
    yMean <- mean(y, na.rm=TRUE);
    segs[idxs,"tcn.mean"] <- rep(yMean, times=length(idxs));

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  verbose && enter(verbose, "Updating DH fields");
  for (kk in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", kk, nbrOfSegments));

    segsKK <- segs[kk,,drop=FALSE];

    # Identify loci in segment
    start <- segsKK$dh.loc.start[1];
    stop <- segsKK$dh.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    # Identify SNPs and non-SNPs
    muN <- data$muN[units];
    isSnp <- is.finite(muN);

    # Identify heterozygous SNPs
    isHet <- (muN == 1/2);
    hets <- whichVector(isSnp & isHet);

    # Update DH mean level
    beta <- data[[by]][units];
    y <- 2*abs(beta[isHet] - 1/2);
    yMean <- mean(y, na.rm=TRUE);
    segs[kk,"dh.mean"] <- yMean;

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Store results
  fitB <- fit;
  fitB$data <- data;
  fitB$output <- segs;

  verbose && exit(verbose);

  fitB;
}) # resampleA()



setMethodS3("resampleB", "PairedPSCBS", function(fit, by=c("betaTN", "betaT"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'by':
  by <- match.arg(by);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Resample (TCN,BAF) signals and reestimate segmentation means");


  fields <- c("CT", "betaN", "betaT", "betaTN", "muN");
  fields <- intersect(names(data), fields);
  verbose && cat(verbose, "Fields to be resampled:");
  verbose && print(verbose, fields);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample TCN segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- fit$data;
  x <- data$x;
  segs <- fit$output;
  nbrOfSegments <- nrow(segs);

  verbose && enter(verbose, "Resampling TCN segments");
  ids <- unique(segs[["tcn.id"]]);
  for (ii in seq(along=ids)) {
    id <- ids[ii];
    verbose && enter(verbose, sprintf("TCN segment #%d of %d", ii, length(ids)));

    # Identify all subsegments for this TCN segment
    idxs <- which(segs[["tcn.id"]] == id);
    n <- length(idxs);

    segsII <- segs[idxs,,drop=FALSE];

    # Identify loci in segment
    start <- segsII$tcn.loc.start[1];
    stop <- segsII$tcn.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    if (length(units) > 1) {
      # Resample
      unitsS <- sample(units, replace=TRUE);
  
      # Resample data
      for (ff in seq(along=fields)) {
        data[[ff]][units] <- data[[ff]][unitsS];
      }
    }

    verbose && exit(verbose);
  } # for (ii ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample DH segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (jj in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("TCN segment #%d of %d", jj, nbrOfSegments));
    segsJJ <- segs[jj,,drop=FALSE];

    # Identify loci in segment
    start <- segsJJ$dh.loc.start[1];
    stop <- segsJJ$dh.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    # Identify SNPs and non-SNPs
    muN <- data$muN[units];
    isSnp <- is.finite(muN);

    # Identify heterozygous SNPs
    isHet <- (muN == 1/2);
    hets <- whichVector(isSnp & isHet);

    units <- units[hets];

    if (length(units) > 1) {
      # Resample
      unitsS <- sample(units, replace=TRUE);
  
      # Resample data
      for (ff in seq(along=fields)) {
        data[[ff]][units] <- data[[ff]][unitsS];
      }
    }

    verbose && exit(verbose);
  } # for (jj ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the segmentation parameter estimates, e.g. mean levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  x <- data$x;

  verbose && enter(verbose, "Updating TCN fields");
  ids <- unique(segs[["tcn.id"]]);
  for (kk in seq(along=ids)) {
    id <- ids[kk];
    verbose && enter(verbose, sprintf("TCN segment #%d of %d", kk, length(ids)));

    # Identify all subsegments for this TCN segment
    idxs <- which(segs[["tcn.id"]] == id);
    n <- length(idxs);

    segsKK <- segs[idxs,,drop=FALSE];

    # Identify loci in segment
    start <- segsKK$tcn.loc.start[1];
    stop <- segsKK$tcn.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    # Update TCN mean level
    y <- data$CT[units];
    yMean <- mean(y, na.rm=TRUE);
    segs[idxs,"tcn.mean"] <- rep(yMean, times=length(idxs));

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  verbose && enter(verbose, "Updating DH fields");
  for (kk in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", kk, nbrOfSegments));

    segsKK <- segs[kk,,drop=FALSE];

    # Identify loci in segment
    start <- segsKK$dh.loc.start[1];
    stop <- segsKK$dh.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    # Identify SNPs and non-SNPs
    muN <- data$muN[units];
    isSnp <- is.finite(muN);

    # Identify heterozygous SNPs
    isHet <- (muN == 1/2);
    hets <- whichVector(isSnp & isHet);

    # Update DH mean level
    beta <- data[[by]][units];
    y <- 2*abs(beta[isHet] - 1/2);
    yMean <- mean(y, na.rm=TRUE);
    segs[kk,"dh.mean"] <- yMean;

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Store results
  fitB <- fit;
  fitB$data <- data;
  fitB$output <- segs;

  verbose && exit(verbose);

  fitB;
}) # resampleB()

##############################################################################
# HISTORY
# 2010-09-16 [HB]
# o Added bootstrap().
# o Added resampleC().
# 2010-09-15 [HB]
# o Added resampleA() and resampleB().  Both are probably incorrect for
#   bootstrapping.
##############################################################################
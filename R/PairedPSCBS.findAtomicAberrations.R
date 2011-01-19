setMethodS3("findAtomicAberrations", "PairedPSCBS", function(this, H=1, alpha=0.02, flavor=c("mean(tcn)", "t(tcn)", "mean(c1,c2)"), verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'H':
  H <- Arguments$getInteger(H, range=c(0,Inf));

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfChromosomes <-nbrOfChromosomes(this);
  if (nbrOfChromosomes > 1) {
    throw("More than one chromosome.");
  }

  nbrOfSegments <- nbrOfSegments(this);

  # Nothing to do?
  if (nbrOfSegments < H+2) {
    res <- list(
      atomicRegions=integer(0),
      atomicIslands=integer(0)
    );
    return(res);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Which dimensions should be tested?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if (flavor == "mean(tcn)") {
    extractSignals <- function(..., na.rm=TRUE) {
      C <- extractLocusLevelTCN(...);  
      # Drop missing values?
      if (na.rm) {
        C <- C[is.finite(C)];
      }
      C;
    } # extractSignals()
    testEquality <- testEqualityTcnByMean;
  } else if (flavor == "t(tcn)") {
    extractSignals <- function(..., na.rm=TRUE) {
      C <- extractLocusLevelTCN(...);  
      # Drop missing values?
      if (na.rm) {
        C <- C[is.finite(C)];
      }
      C;
    } # extractSignals()
    testEquality <- testEqualityTcnByT;
  } else if (flavor == "mean(c1,c2)") {
    extractSignals <- function(..., na.rm=TRUE) {
      data <- extractLocusLevelC1C2(...);  
      # Drop missing values?
      if (na.rm) {
        ok <- is.finite(data$C1) & is.finite(data$C2);
        data <- data[ok,,drop=FALSE];
      }
      data <- as.matrix(data[,c("C1","C2")]);
      data;
    } # extractSignals()
    testEquality <- testEqualityC1C2ByMean;
  }

  verbose && enter(verbose, "Call equivalent copy-number states by pruning");

  # Initial set of atomic regions
  atomicRegions <- NULL;

  nbrOfAberrations <- (nbrOfSegments-H);
  for (rr in 2:nbrOfAberrations) {
    verbose && enter(verbose, sprintf("Aberration #%d of %d", rr, nbrOfAberrations));
    verbose && printf(verbose, "alpha=%f\n", alpha);

    # The two flanking regions
    idxL <- rr-1;
    idxR <- rr+H;
    fitL <- extractByRegion(this, region=idxL);
    fitR <- extractByRegion(this, region=idxR);

    # Extract their data
    verbose && enter(verbose, "Extracting signals");
    dataL <- extractSignals(fitL);
    dataR <- extractSignals(fitR);
    verbose && str(verbose, dataL);
    verbose && str(verbose, dataR);
    verbose && exit(verbose);

    # Test if they are equal
    isEqual <- testEquality(dataL, dataR, alpha=alpha);
    verbose && print(verbose, isEqual);
    fit <- attr(isEqual, "fit");
    # Drop attributes
    isEqual <- as.logical(isEqual);

    verbose && printf(verbose, "t=%.3f (p=%g), (alpha=%g) (L==R)=%s\n", 
                       fit$statistic, fit$p.value, alpha, isEqual);
    rm(dataL, dataR, fit); # Not needed anymore

    # If the two flanking regions are equal, then we have 
    # found an atomic region.
    if (isTRUE(isEqual)) {
      atomicRegions <- c(atomicRegions, rr);
      verbose && print(verbose, atomicRegions);
    }

    verbose && exit(verbose);
  } # for (rr ...)

  # Table of atomic regions of length K found
  res <- data.frame(
    leftRegion  = atomicRegions-1L,
    rightRegion = atomicRegions+(H-1L)+1L,
    firstRegion = atomicRegions,
    lastRegion  = atomicRegions+(H-1L)
#    start       = start[atomicRegions],
#    stop        = stop[atomicRegions+(H-1L)]
  );

  # Atomic islands = atomic regions that are not next 
  # to another atomic region
  dups <- which(diff(atomicRegions) == 1);
  if (length(dups) > 0) {
    dups <- c(dups, dups+1L);
    atomicIslands <- atomicRegions[-dups];
  } else {
    atomicIslands <- atomicRegions;
  }

  res <- list(
    H=H,
    atomicRegions=atomicRegions,
    atomicIslands=atomicIslands,
    ambigousRegions=setdiff(atomicRegions, atomicIslands),
    res=res
  );

  verbose && exit(verbose);

  res;
}, protected=TRUE) # findAtomicAberrations()


############################################################################
# HISTORY:
# 2011-01-12
# o Created findAtomicAberrations() for PairedPSCBS from ditto for
#   CopyNumberRegions currently in aroma.cn.
# o Created.
############################################################################

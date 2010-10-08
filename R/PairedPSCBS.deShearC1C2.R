###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod deShearC1C2
#
# @title "Correct for shearing in (C1,C2) space based on region-based PSCN estimates"
#
# \description{
#  @get "title" as given by the PSCBS segmentation method.
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A PairedPSCBS fit object as returned by 
#     @see "psCBS::segmentByPairedPSCBS".}
#   \item{adjust}{A @numeric adjusting the bandwidth of the empirical
#     density estimator of line (changepoint) directions.}
#   \item{weightFlavor}{A @character string specifying how weights are
#     generated for lines (changepoints).}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a PairedPSCBS fit object.
# }
#
# @examples "../incl/deShearC1C2.PairedPSCBS.Rex"
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("deShearC1C2_v0", "PairedPSCBS", function(fit, adjust=0.5, tol=0.02, flavor=c("decreasing", "all"), weightFlavor=c("min", "sum"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0,Inf));

  flavor <- match.arg(flavor);

  # Argument 'weightFlavor':
  weightFlavor <- match.arg(weightFlavor);
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Correct for shearing in (C1,C2) space at the region level ");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  stopifnot(!is.null(data));

  segs <- as.data.frame(fit);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # (C1,C2,...)
  X <- extractC1C2(fit);

  # Number of TCN and DH data points
  counts <- X[,3:4, drop=FALSE];

  # Region weights from DH counts
  w <- counts[,2];
  w <- sqrt(w);
  w <- w / sum(w, na.rm=TRUE);

  ## Change-point weights
  if (weightFlavor == "min") {
    ## smallest of the two flanking (DH) counts 
    cpw <- cbind(w[1:(length(w)-1)], w[2:length(w)]);
    cpw <- rowMins(cpw, na.rm=TRUE);
    cpw[is.infinite(cpw)] <- NA;
    cpw <- sqrt(cpw);
  } else if (weightFlavor == "sum") {
    ## sum of region weights
    cpw <- w[1:(length(w)-1)] + w[2:length(w)];
  }
  
  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];
  dC1C2 <- colDiffs(C1C2);
  alpha <- atan(dC1C2[,2]/dC1C2[,1]);
  radius <- sqrt(dC1C2[,2]^2 + dC1C2[,1]^2);

  verbose && enter(verbose, "Represent change points as angles");
  verbose && str(verbose, alpha);

  # Keep only finite data points
  ok <- (is.finite(alpha) & is.finite(cpw));
  alphaT <- alpha[ok];
  verbose && cat(verbose, "Finite slopes in (C1,C2):");
  cpwT <- cpw[ok];
  cpwT <- cpwT / sum(cpwT);
  verbose && print(verbose, cbind(alphaT=alphaT, wT=cpwT));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Adjust modes in {alpha} to their expect locations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Shift, modulo pi/8:th, and unshift");
  ## The signal is pi-periodic.
  ## we are looking at from -pi/2 to pi/2.
  ## we expect a peak near -pi/2 (or pi/2...)
  ## in order to estimate it correctly, transform the signal so that it is
  ## in -pi/2-pi/8, pi/2-pi/8
  ## /PN 2010-09-22
  lag <- pi/8;
  aa <- alphaT;
  ww <- which(alphaT > pi/2-lag) ## half way to both expected peaks
  aa[ww] <- aa[ww]-pi;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Find modes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Finding modes (peaks & valleys)");
  rg <- range(aa);
  fp <- findPeaksAndValleys(aa, weights=cpwT, from=rg[1], to=rg[2], adjust=adjust, tol=tol, ...);
  verbose && cat(verbose, "Peaks and valleys:");
  verbose && print(verbose, fp);

  verbose && cat(verbose, "Peaks:");
  type <- NULL; rm(type); # To please R CMD check
  pfp <- subset(fp, type == "peak");
  verbose && print(verbose, pfp);
  nbrOfPeaks <- nrow(pfp);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Call modes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calling modes");
  expected <- c(-1/2,-1/4,0,+1/4,+1/2)*pi;
  pfp <- callPeaks(pfp, expected=expected, flavor=flavor, verbose=verbose);
  verbose && cat(verbose, "Calls:");
  verbose && print(verbose, pfp);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Estimate (C1,C2) shear model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Estimate (C1,C2) shear model");
  ## Use -pi/2 and 0 to correct for shearing

  # (a) Vertical shear (based on horizontal information)
  pfpT <- subset(pfp, call == -pi/2);
  # Sanity checks
  stopifnot(nrow(pfpT) == 1);
  # Shear parameter
  tX <- pfpT$x;
#  residual <- (pfpT$x - -pi/2);
#  tX <- pi - residual;
  # Sanity checks
  stopifnot(is.finite(tX));

  # (b) Horizontal shear (based on vertical information)
  pfpT <- subset(pfp, call == 0);
  # Sanity checks
  stopifnot(nrow(pfpT) == 1);
  # Shear parameter
  tY <- pi/2 - pfpT$x;
#  residual <- (pfpT$x - 0);
#  tY <- pi - residual;
  # Sanity checks
  stopifnot(is.finite(tY));

  # Create backtransform function
  H <- function(xy) {
    sx <- tan(tX+pi/2);
    sy <- tan(tY+pi/2);
#    sx <- tan(tX);
#    sy <- tan(tY);
    cbind(xy[, 1]+xy[, 2]*sx, xy[, 1]*sy + xy[, 2]);
  } # H()

  verbose && cat(verbose, "Model fit:");
  modelFit <- list(H=H, parameters=c(tX=tX, tY=tY), debug=list(pfp=pfp));
  verbose && str(verbose, modelFit);

  # Not needed anymore
  rm(pfpT);
  verbose && exit(verbose);



  verbose && enter(verbose, "Remove shearing from the (C1,C2) space by independent linear adjustment in x and y");
  ## Backtransform
  C1C2o <- H(C1C2);
  # Sanity checks
  stopifnot(dim(C1C2o) == dim(C1C2));
  verbose && exit(verbose);

  verbose && enter(verbose, "Transform orthogonalized (C1,C2) to (TCN,DH)");
  # (C1,C2) -> (TCN,DH)
  gamma <- rowSums(C1C2o, na.rm=TRUE);
  dh <- 2*(C1C2o[,2]/gamma - 1/2);
  verbose && exit(verbose);

  # Update segmentation means
  segs[,"tcn.mean"] <- gamma;
  segs[,"dh.mean"] <- dh;

  # Update data [TO DO]


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fitO <- fit;
  fitO$output <- segs;
  fitO$deshearModelFit <- modelFit;

  verbose && exit(verbose);

  fitO;
}) # deShearC1C2_v0()






setMethodS3("fitC1C2ShearModel", "PairedPSCBS", function(fit, adjust=0.5, tol=0.02, flavor=c("decreasing", "all"), weightFlavor=c("min", "sum"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0,Inf));

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'weightFlavor':
  weightFlavor <- match.arg(weightFlavor);
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting (C1,C2) shear model");

  verbose && enter(verbose, "Extracting (C1,C2) => (radius,alpha)");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  stopifnot(!is.null(data));

  segs <- as.data.frame(fit);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # (C1,C2,...)
  X <- extractC1C2(fit);

  # Number of TCN and DH data points
  counts <- X[,3:4, drop=FALSE];

  # Region weights from DH counts
  w <- counts[,2];
  w <- sqrt(w);
  w <- w / sum(w, na.rm=TRUE);

  ## Change-point weights
  if (weightFlavor == "min") {
    ## smallest of the two flanking (DH) counts 
    cpw <- cbind(w[1:(length(w)-1)], w[2:length(w)]);
    cpw <- rowMins(cpw, na.rm=TRUE);
    cpw[is.infinite(cpw)] <- NA;
    cpw <- sqrt(cpw);
  } else if (weightFlavor == "sum") {
    ## sum of region weights
    cpw <- w[1:(length(w)-1)] + w[2:length(w)];
  }
  
  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];
  dC1C2 <- colDiffs(C1C2);

  # (radius,alpha)  # 'radius' is not really used
  radius <- sqrt(dC1C2[,2]^2 + dC1C2[,1]^2);
  alpha <- atan(dC1C2[,2]/dC1C2[,1]);

  verbose && cat(verbose, "(radius,alpha):");
  verbose && str(verbose, radius);
  verbose && str(verbose, alpha);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Dropping non-finite signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Dropping non-finite signals");
  # Keep only finite data points
  ok <- (is.finite(radius) & is.finite(alpha) & is.finite(cpw));
  radiusT <- radius[ok];
  alphaT <- alpha[ok];
  verbose && cat(verbose, "Finite slopes in (C1,C2):");
  cpwT <- cpw[ok];
  cpwT <- cpwT / sum(cpwT);
  verbose && cat(verbose, "(radius,alpha,weights):");
  verbose && print(verbose, cbind(radiusT=radiusT, alphaT=alphaT, wT=cpwT));
  rm(ok);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Adjust modes in {alpha} to their expect locations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Shift, modulo pi/8:th, and unshift");
  ## The signal is pi-periodic.
  ## we are looking at from -pi/2 to pi/2.
  ## we expect a peak near -pi/2 (or pi/2...)
  ## in order to estimate it correctly, transform the signal so that it is
  ## in -pi/2-pi/8, pi/2-pi/8
  ## /PN 2010-09-22
  lag <- pi/8;
  aa <- alphaT;
  ww <- which(alphaT > pi/2-lag) ## half way to both expected peaks
  aa[ww] <- aa[ww]-pi;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Find modes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Finding modes (peaks & valleys)");
  rg <- c(-pi/2,pi/2)-pi/8;
  fp <- findPeaksAndValleys(aa, weights=cpwT, from=rg[1], to=rg[2], adjust=adjust, tol=tol, ...);
  verbose && cat(verbose, "Peaks and valleys:");
  verbose && print(verbose, fp);

  verbose && cat(verbose, "Peaks:");
  type <- NULL; rm(type); # To please R CMD check
  pfp <- subset(fp, type == "peak");
  verbose && print(verbose, pfp);
  nbrOfPeaks <- nrow(pfp);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Call modes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calling modes");
  expected <- c(-1/2,-1/4,0,+1/4,+1/2)*pi;
  verbose && cat(verbose, "Expected locations of peaks:");
  verbose && print(verbose, expected);
  pfp <- callPeaks(pfp, expected=expected, flavor=flavor, verbose=verbose);
  verbose && cat(verbose, "Calls:");
  verbose && print(verbose, pfp);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Estimate (C1,C2) shear model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Estimate (C1,C2) shear model");
  ## Use -pi/2 and 0 to correct for shearing

  # (a) Vertical shear (based on horizontal information)
  pfpT <- subset(pfp, call == -pi/2);
  # Sanity checks
  stopifnot(nrow(pfpT) == 1);
  # Shear parameter
  tX <- pfpT$x;
#  residual <- (pfpT$x - -pi/2);
#  tX <- pi - residual;
  # Sanity checks
  stopifnot(is.finite(tX));

  # (b) Horizontal shear (based on vertical information)
  pfpT <- subset(pfp, call == 0);
  # Sanity checks
  stopifnot(nrow(pfpT) == 1);
  # Shear parameter
  tY <- pi/2 - pfpT$x;
#  residual <- (pfpT$x - 0);
#  tY <- pi - residual;
  # Sanity checks
  stopifnot(is.finite(tY));

  # Create backtransform function
  H <- function(xy) {
    sx <- tan(tX+pi/2);
    sy <- tan(tY+pi/2);
#    sx <- tan(tX);
#    sy <- tan(tY);
    cbind(xy[, 1]+xy[, 2]*sx, xy[, 1]*sy + xy[, 2]);
  } # H()

  verbose && cat(verbose, "Model fit:");
  modelFit <- list(H=H, parameters=c(tX=tX, tY=tY), debug=list(pfp=pfp));
  verbose && str(verbose, modelFit);

  # Not needed anymore
  rm(pfpT);
  verbose && exit(verbose);

  verbose && exit(verbose);

  modelFit;
}) # fitC1C2ShearModel()



setMethodS3("deShearC1C2", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Correct for shearing in (C1,C2) space at the region level ");

  verbose && enter(verbose, "Fitting (C1,C2) shear model");
  modelFit <- fitC1C2ShearModel(fit, ..., verbose=verbose);
  verbose && print(verbose, modelFit);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  stopifnot(!is.null(data));

  segs <- as.data.frame(fit);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # (C1,C2,...)
  X <- extractC1C2(fit);

  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];


  verbose && enter(verbose, "Remove shearing from the (C1,C2) space by independent linear adjustment in x and y");
  ## Backtransform
  H <- modelFit$H;

  C1C2o <- H(C1C2);
  # Sanity checks
  stopifnot(dim(C1C2o) == dim(C1C2));
  verbose && exit(verbose);

  verbose && enter(verbose, "Transform orthogonalized (C1,C2) to (TCN,DH)");
  # (C1,C2) -> (TCN,DH)
  gamma <- rowSums(C1C2o, na.rm=TRUE);
  dh <- 2*(C1C2o[,2]/gamma - 1/2);
  verbose && exit(verbose);

  # Update segmentation means
  segs[,"tcn.mean"] <- gamma;
  segs[,"dh.mean"] <- dh;

  # Update data [TO DO]


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fitO <- fit;
  fitO$output <- segs;
  fitO$deshearModelFit <- modelFit;

  verbose && exit(verbose);

  fitO;
}) # deShearC1C2()


##############################################################################
# HISTORY
# 2010-10-08 [HB]
# o Now deShearC1C2() uses fitC1C2ShearModel().
# o Added fitC1C2ShearModel().
# o Now deShearC1C2() returns the 'modelFit'.
# o Now deShearC1C2() calls peaks using callPeaks() for PeaksAndValleys.
# 2010-09-26 [HB]
# o Added argument 'adjust' to deShearC1C2() with new default.
# o Added sanity checks to deShearC1C2().
# o Now normalizeBAFsByRegions() for PairedPSCBS handles multiple chromosomes.
# 2010-09-22 [PN]
# o Added deShearC1C2() for PairedPSCBS.
# 2010-09-19 [HB+PN]
# o Added orthogonalizeC1C2() for PairedPSCBS.
# 2010-09-15 [HB]
# o Added Rdocs for callCopyNeutralRegions().
# 2010-09-09 [HB]
# o Added callCopyNeutralRegions() for PairedPSCBS.
# 2010-09-08 [HB]
# o Added subsetBySegments() for PairedPSCBS.
# o Added Rdocs with an example.
# o Added normalizeBAFsByRegions() for PairedPCSBS.
# o Created.
##############################################################################

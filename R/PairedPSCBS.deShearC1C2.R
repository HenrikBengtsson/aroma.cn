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
setMethodS3("deShearC1C2", "PairedPSCBS", function(fit, ..., dirs=c("|-", "-", "|", "X", "|,-", "-,|", "|-,X", "|,-,X", "-,|,X"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'dirs':
  dirs <- match.arg(dirs);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Correct for shearing in (C1,C2) space at the region level ");
  verbose && cat(verbose, "Directions: ", paste(dirs, collapse=","));

  # First one dimension, then the other?
  dirs <- unlist(strsplit(dirs, split=",", fixed=TRUE));
  if (length(dirs) > 1) {
    verbose && enter(verbose, "Deshearing direction by direction");
    fitT <- fit;
    for (kk in seq(along=dirs)) {
      dir <- dirs[kk];
      verbose && enter(verbose, sprintf("Direction #%d ('%s') of %d", kk, dir, length(dirs)));
      fitT <- deShearC1C2(fitT, dirs=dir, ..., verbose=verbose);
      verbose && exit(verbose);
    }
    verbose && exit(verbose);

    verbose && exit(verbose);
    return(fitT);
  }

  verbose && enter(verbose, "Fitting (C1,C2) shear model");
  modelFit <- fitDeltaC1C2ShearModel(fit, ..., verbose=verbose);
  verbose && print(verbose, modelFit);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  if (dirs == "|-") {
    H <- modelFit$H;
  } else if (dirs == "-") {
    H <- modelFit$Hx;
  } else if (dirs == "|") {
    H <- modelFit$Hy;
  } else if (dirs == "X") {
    H <- modelFit$Hd;
  }

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
  fitO$modelFit <- modelFit;

  verbose && exit(verbose);

  fitO;
}) # deShearC1C2()



setMethodS3("translateC1C2", "PairedPSCBS", function(fit, dC1=0, dC2=0, sC1=1, sC2=1, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dC1 <- Arguments$getNumeric(dC1);
  dC2 <- Arguments$getNumeric(dC2);
  sC1 <- Arguments$getNumeric(sC1);
  sC2 <- Arguments$getNumeric(sC2);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Scaling and translating by (C1,C2)");
  verbose && cat(verbose, "sC1: ", sC1);
  verbose && cat(verbose, "sC2: ", sC2);
  verbose && cat(verbose, "dC1: ", dC1);
  verbose && cat(verbose, "dC2: ", dC2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- as.data.frame(fit);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # (C1,C2,...)
  X <- extractC1C2(fit);

  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];

  C1C2[,1] <- dC1 + sC1*C1C2[,1];
  C1C2[,2] <- dC2 + sC2*C1C2[,2];

  verbose && enter(verbose, "(C1,C2) to (TCN,DH)");
  # (C1,C2) -> (TCN,DH)
  gamma <- rowSums(C1C2, na.rm=TRUE);
  dh <- 2*(C1C2[,2]/gamma - 1/2);
  verbose && exit(verbose);

  # Update segmentation means
  segs[,"tcn.mean"] <- gamma;
  segs[,"dh.mean"] <- dh;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fitO <- fit;
  fitO$output <- segs;

  verbose && exit(verbose);

  fitO;
})



setMethodS3("transformC1C2", "PairedPSCBS", function(fit, fcn, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fcn':
  stopifnot(is.function(fcn));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Transform (C1,C2) by a function");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- as.data.frame(fit);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # (C1,C2,...)
  X <- extractC1C2(fit);

  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];

  C1C2 <- fcn(C1C2, ...);

  verbose && enter(verbose, "(C1,C2) to (TCN,DH)");
  # (C1,C2) -> (TCN,DH)
  gamma <- rowSums(C1C2, na.rm=TRUE);
  dh <- 2*(C1C2[,2]/gamma - 1/2);
  verbose && exit(verbose);

  # Update segmentation means
  segs[,"tcn.mean"] <- gamma;
  segs[,"dh.mean"] <- dh;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fitO <- fit;
  fitO$output <- segs;

  verbose && exit(verbose);

  fitO;
})



setMethodS3("fitDeltaC1C2ShearModel", "PairedPSCBS", function(fit, adjust=0.5, tol=0.02, flavor=c("decreasing", "all"), weightFlavor=c("min", "sum"), ..., verbose=FALSE) {
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


  verbose && enter(verbose, "Fitting (C1,C2) shear model by (dC1,dC2)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- as.data.frame(fit);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # (C1,C2,...)
  X <- extractC1C2(fit);

  # Number of TCN and DH data points
  counts <- X[,3:4, drop=FALSE];

  # (C1,C2)
  X <- X[,1:2,drop=FALSE];
  X <- as.matrix(X);

  # Region weights from DH counts
  w <- counts[,2];
  w <- sqrt(w);
  rm(counts);

  ## Change-point weights
  if (weightFlavor == "min") {
    ## smallest of the two flanking (DH) counts 
    w <- cbind(w[1:(length(w)-1)], w[2:length(w)]);
    w <- rowMins(w, na.rm=TRUE);
    w[is.infinite(w)] <- NA;
    dw <- sqrt(w);
  } else if (weightFlavor == "sum") {
    ## sum of region weights
    dw <- w[1:(length(w)-1)] + w[2:length(w)];
  }


  modelFit <- fitDeltaXYShearModel(X, weights=dw, adjust=adjust, ..., verbose=less(verbose, 1));

  verbose && cat(verbose, "Model fit:");
  verbose && str(verbose, modelFit);

  verbose && exit(verbose);

  modelFit;
}) # fitDeltaC1C2ShearModel()



setMethodS3("fitDeltaXYShearModel", "matrix", function(X, weights=NULL, adjust=0.5, tol=0.02, flavor=c("decreasing", "all"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'X':
  dim <- dim(X);
  if (dim[2] != 2) {
    throw("Argument 'X' is not a Jx2 matrix: ", paste(dim, collapse="x"));
  }
  J <- dim[1];

  # Argument 'weights':
  length2 <- rep(J-1L, times=2);
  if (!is.null(weights)) {
    weights <- Arguments$getDoubles(weights, range=c(0,Inf), 
                                       length=length2, disallow=NULL);
  }

  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0,Inf));

  # Argument 'flavor':
  flavor <- match.arg(flavor);
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting (x,y) shear model");

  verbose && enter(verbose, "Extracting (radius,alpha)");
  verbose && cat(verbose, "Number of (x,y) data points: ", J);
  verbose && cat(verbose, "Number of (dx,dy) data points: ", J-1L);

  # (dX,dY)
  dXY <- colDiffs(X);

  # (radius,alpha)  # 'radius' is not really used
  radius <- sqrt(dXY[,2]^2 + dXY[,1]^2);
  alpha <- atan(dXY[,2]/dXY[,1]);

  verbose && cat(verbose, "(radius,alpha):");
  verbose && str(verbose, radius);
  verbose && str(verbose, alpha);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Dropping non-finite signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Dropping non-finite signals");
  # Keep only finite data points
  ok <- (is.finite(radius) & is.finite(alpha) & is.finite(weights));
  radius <- radius[ok];
  alpha <- alpha[ok];
  verbose && cat(verbose, "Finite *transitions* in (C1,C2):");
  weights <- weights[ok];

  # Standardize weights to [0,1]
  weights <- weights / sum(weights, na.rm=TRUE);

  verbose && cat(verbose, "(radius,alpha,weights):");
  verbose && print(verbose, cbind(radius=radius, alpha=alpha, weights=weights));
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
  alphaT <- alpha;
  ww <- which(alphaT > pi/2-lag) ## half way to both expected peaks
  alphaT[ww] <- alphaT[ww]-pi;
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Find modes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Finding modes (peaks & valleys)");
  rg <- c(-pi/2,pi/2)-pi/8;
  d <- density(alphaT, weights=weights, from=rg[1], to=rg[2], adjust=adjust);

  fp <- findPeaksAndValleys(d, tol=tol, ...);
  # Not needed anymore
  rm(alphaT, rg);

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
  expected <- c("-"=-1/2, "\\"=-1/4, "|"=0, "/"=+1/4, "-"=+1/2)*pi;

  verbose && cat(verbose, "Expected locations of peaks:");
  verbose && print(verbose, expected);

  pfp <- callPeaks(pfp, expected=expected, flavor=flavor, verbose=verbose);

  verbose && cat(verbose, "Calls:");
  verbose && print(verbose, pfp);
  verbose && exit(verbose);

  if (TRUE) {
    devSet(sprintf("d,%s", digest(d)));
    plot(d, lwd=2);
    abline(v=expected);
    text(x=expected, y=par("usr")[4], names(expected), adj=c(0.5,-0.5), cex=1, xpd=TRUE);
    idxs <- match(pfp$call, expected);
    text(x=pfp$x, y=pfp$density, names(expected)[idxs], adj=c(0.5, -0.5), cex=1, col="blue");
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Estimate (x,y) shear model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Estimate (x,y) shear model");
  ## Use -pi/2 and 0 to correct for shearing

  # (a) Vertical shear (based on horizontal information)
  pfpT <- subset(pfp, call == -pi/2);
  verbose && cat(verbose, "Horizontal peak:");
  verbose && print(verbose, pfpT);
  # Sanity checks
  if (nrow(pfpT) < 1) {
    throw(sprintf("Cannot fit vertical shear parameter. No lines where called to horizontal (angle=-pi/2=%f).", -pi/2));
  }
  # Shear parameter
  tX <- pfpT$x;
#  residual <- (pfpT$x - -pi/2);
#  tX <- pi - residual;
  # Sanity checks
  stopifnot(is.finite(tX));
  # Shear parameter
  sx <- tan(tX+pi/2);
  stopifnot(is.finite(sx));

  # (b) Horizontal shear (based on vertical information)
  pfpT <- subset(pfp, call == 0);
  verbose && cat(verbose, "Vertical peak:");
  verbose && print(verbose, pfpT);
  # Sanity checks
  if (nrow(pfpT) < 1) {
    throw("Cannot fit horizontal shear parameter. No lines where called to vertical (angle=0).");
  }
  # Shear parameter
  tY <- pi/2 - pfpT$x;
#  residual <- (pfpT$x - 0);
#  tY <- pi - residual;
  # Sanity checks
  stopifnot(is.finite(tY));
  # Shear parameter
  sy <- tan(tY+pi/2);
  stopifnot(is.finite(sy));

  # (c) Vertical scale (based on diagonal information)
  pfpT <- subset(pfp, call %in% c(-pi/4, +pi/4));
  verbose && cat(verbose, "Diagonal peak:");
  verbose && print(verbose, pfpT);
  # Sanity checks
  if (nrow(pfpT) < 1) {
    msg <- "Cannot fit diagonal shear scale parameter. No lines where called to diagonal (angle=-pi/4 or +pi/4).";
    warning(msg);

    scale <- scaleY <- as.double(NA);
    Hd <- NULL;
  } else {
    # Shear parameter
    tX <- pfpT$x;
    scaleY <- pfpT$call/pfpT$x;
    verbose && cat(verbose, "scaleY:");
    verbose && print(verbose, scaleY);
    scaleY <- weighted.mean(scaleY, w=pfpT$density);
    verbose && print(verbose, scaleY);
  
    # Preserve (dx^2 + dy^2) before and after
    scale <- 1/sqrt(1+sx^2);

    Hd <- function(xy) {
      y <- xy[,2];
      mu0 <- min(y, na.rm=TRUE);
      y <- scaleY * y;
      mu1 <- min(y, na.rm=TRUE);
      dmu <- mu1 - mu0;
      # Preserve y offset
      y <- y - dmu;
      xy[,2] <- y;
      xy;
    } # Hd()
  }
  
  # Create backtransform function
  H <- function(xy) {
    cbind(xy[, 1]+xy[, 2]*sx, xy[, 1]*sy + xy[, 2]);
  } # H()


  Hx <- function(xy) {
    cbind(xy[, 1]+xy[, 2]*sx, xy[, 2]);
  } # Hx()

  Hy <- function(xy) {
    cbind(xy[, 1], xy[, 1]*sy + xy[, 2]);
  } # Hy()


  # Alternative?!?
  T <- matrix(c(1,sx, sy,1), nrow=2, ncol=2, byrow=FALSE);
  H2 <- function(xy) {
    xy <- as.matrix(xy);
    xy %*% t(T);
  } # H2()

  verbose && cat(verbose, "Model fit:");
  modelFit <- list(H=H, Hx=Hx, Hy=Hy, Hd=Hd, parameters=c(sx=sx, sy=sy, scaleY=scaleY, scale=scale, tX=tX, tY=tY), debug=list(pfp=pfp));
  verbose && str(verbose, modelFit);

  # Sanity checks
  range <- c(-3,3);
  sx <- Arguments$getDouble(sx, range=range);
  sy <- Arguments$getDouble(sy, range=range);

  # Not needed anymore
  rm(pfpT);
  verbose && exit(verbose);

  verbose && exit(verbose);

  modelFit;
}) # fitDeltaXYShearModel()


setMethodS3("estimateC2Bias", "PairedPSCBS", function(fit, ...) {
  # Identify region in allelic balance
  segs <- as.data.frame(fit);
  ab.call <- segs$ab.call;
  if (is.null(ab.call)) {
    throw("Allelic balance has not been called.");
  }
  idxs <- which(segs$ab.call);
  segs <- segs[idxs,,drop=FALSE];

  # Extract (TCN,DH)
  gamma <- segs[, "tcn.mean"];
  rho <- segs[, "dh.mean"];

  # Calculate (C1,C2)
  C1 <- 1/2 * (1 - rho) * gamma;
  C2 <- gamma - C1;

  # Calculate bias in C2 
  dC2 <- C2 - C1;

  # Calculate weighted average of all C2 biases
  w <- segs[,"dh.num.mark"];
  w <- w / sum(w, na.rm=TRUE);
  dC2 <- weightedMedian(dC2, w=w);

  dC2;
}) # estimateC2Bias()



##############################################################################
# HISTORY
# 2010-10-20 [HB]
# o Now fitDeltaXYShearModel() uses both -pi/4 and +pi/4 to estimate
#   the diagonal parameters.
# 2010-10-08 [HB]
# o Added estimateC2Bias().
# o Added fitDeltaXYShearModel().
# o Now deShearC1C2() uses fitDeltaC1C2ShearModel().
# o Added fitDeltaC1C2ShearModel().
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

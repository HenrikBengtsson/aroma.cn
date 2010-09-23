
  C1C2toAB <- function(xy, ...) {
    dxy <- colDiffs(xy[,1:2]);
#print(dxy);
    # (l,b) = Length and slope of line
    l <- sqrt(dxy[,1]^2 + dxy[,2]^2);
    k <- (dxy[,2]/dxy[,1]);
    b <- atan(1/k);
#print(list(l=l, b=b, k=k, dxy=dxy));


    # Regions weights
    counts <- xy[,"dh.num.mark", drop=TRUE];
    w <- sqrt(counts);
    w <- w / sum(w, na.rm=TRUE);

    # Change-point weights
    w <- w[1:(length(w)-1)] + w[2:length(w)];

    res <- as.matrix(cbind(l=l, b=b, w=w, k=k));
#    attr(res, "dxy") <- dxy;
    res;
  } # C1C2toAB()

  ABtoC1C2 <- function(AB, C1C2, ...) {
    C1C2 <- C1C2[,c("C1","C2"), drop=FALSE];
    C1C2 <- as.matrix(C1C2);
    dxy <- colDiffs(C1C2);

    #  b = atan(dC1C2[,1]/dC1C2[,2]);
    #  tan(b) = dC1/dC2;
    #  dC2*tan(b) = dC1;
    l <- AB[,"l"];
    b <- AB[,"b"];
    k <- AB[,"k"];
    dxy2 <- sign(dxy);
    k2 <- tan(b);
    dxy2[,2] <- dxy2[,1]/k2;
    l2 <- sqrt(dxy2[,1]^2+dxy2[,2]^2);
    s <- l/l2;
    dxy2 <- s*dxy2;

    for (cc in 1:2) {
      delta <- dxy2[,cc];
      idxs <- whichVector(is.finite(delta));
      C1C2[idxs+1,cc] <- C1C2[idxs,cc] + delta[idxs];
    }

    C1C2;
  } # ABtoC1C2()

###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod normalizeBAFsByRegions
#
# @title "Normalizes allele B fractions (BAFs) based on region-based PSCN estimates"
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
#   \item{by}{A @character string specifying if the normalization function
#     should be estimated based on TumorBoost normalized or non-normalized
#     tumor allele B fractions (BAFs).}
#   \item{...}{Additional arguments passed
#     @see "aroma.cn::normalizeMirroredBAFsByRegions".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a PairedPSCBS fit object where the region-level 
#   decrease-in-heterozygosity (DH) means have been normalized,
#   as well as the locus-specific tumor allele B fractions.
# }
#
# \details{
#   Note that his normalization method depends on the segmentation
#   results. Hence, it recommended \emph{not} to resegment the
#   normalized signals returned by this, because such a segmentation
#   will be highly dependent on the initial segmentation round.
# }
#
# @examples "../incl/normalizeBAFsByRegions.PairedPSCBS.Rex"
#
# @author
#
# \seealso{
#   Internally @see "aroma.cn::normalizeMirroredBAFsByRegions" is used.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("normalizeBAFsByRegions", "PairedPSCBS", function(fit, by=c("betaTN", "betaT"), ..., verbose=FALSE) {
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


  verbose && enter(verbose, "Normalizes region-level mirrored allele B fractions (mBAFs)");

  data <- fit$data;
  stopifnot(!is.null(data));

  segs <- fit$output;
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  x <- data$x;
  betaT <- data$betaT;
  betaTN <- data$betaTN;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate region-level mBAFs for homozygous SNPs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calculating region-level mBAFs for homozygous SNPs");

  # Calculate mBAFs for all loci
  beta <- data[[by]];
  rho <- 2*abs(beta - 1/2);
  rm(beta);

  # Identify homozygous SNPs
  muN <- data$muN;
  isHom <- (muN == 0 | muN == 1);
  rm(muN);

  # Allocate
  naValue <- as.double(NA);
  X <- matrix(naValue, nrow=nbrOfSegments, ncol=3);
  for (kk in seq(length=nbrOfSegments)) {
    xRange <- as.numeric(segs[kk,c("dh.loc.start", "dh.loc.end")]);
    tcn <- segs[kk,"tcn.mean"];
    dh <- segs[kk,"dh.mean"];
    # Identify all homozygous SNPs in the region
    idxs <- whichVector(xRange[1] <= x & x <= xRange[2] & isHom);
    mBAFhom <- mean(rho[idxs], na.rm=TRUE);
    X[kk,] <- c(dh, mBAFhom, tcn);
  } # for (kk ...)

  # Not needed anymore
  rm(rho, isHom);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Normalize region-level DHs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Normalizing region-level mBAFs");
  XN <- normalizeMirroredBAFsByRegions(data=X, ..., verbose=verbose);

  # Update DH segmentation means
  rhoN <- XN[,1,drop=TRUE];
  segs[,"dh.mean"] <- rhoN;
  rm(rhoN);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Normalize locus-level data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Normalizing locus-level data accordingly");
  modelFit <- attr(XN, "modelFit");
  scale <- modelFit$scale;
  rm(modelFit, XN);

  # Expand region-level scale factors to locus-level scale factors
  naValue <- as.double(NA);
  scales <- rep(naValue, times=length(data$betaT));
  for (kk in seq(length=nbrOfSegments)) {
    xRange <- as.numeric(segs[kk,c("dh.loc.start", "dh.loc.end")]);
    # Identify all SNPs in the region
    idxs <- whichVector(xRange[1] <= x & x <= xRange[2]);
    scales[idxs] <- scale[kk];
  } # for (kk ...)

  # Update tumor allele B fractions
  for (ff in c("betaT", "betaTN")) {
    beta <- data[[ff]];
    beta <- beta - 1/2;
    beta <- scales * beta;
    beta <- beta + 1/2;
    data[[ff]] <- beta;
  }
  rm(beta, scales);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fitN <- fit;
  fitN$data <- data;
  fitN$output <- segs;

  verbose && exit(verbose);

  fitN;
}) # normalizeBAFsByRegions()





###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod orthogonalizeC1C2
#
# @title "Orthogonalizes (C1,C2) based on region-based PSCN estimates"
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
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a PairedPSCBS fit object.
# }
#
# @examples "../incl/orthogonalizeC1C2.PairedPSCBS.Rex"
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("orthogonalizeC1C2", "PairedPSCBS", function(fit, ..., debugPlot=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Orthogonalizes region-level minor and major CNs (C1,C2)");

  data <- fit$data;
  stopifnot(!is.null(data));

  segs <- as.data.frame(fit);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (C1,C2,...)
  X <- extractMinorMajorCNs(fit);

  # Number of TCN and DH data points
  counts <- X[,3:4, drop=FALSE];

  # Regions weights
  w <- counts[,2];
  w <- sqrt(w);
  w <- w / sum(w, na.rm=TRUE);

  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];

  if (debugPlot) {
    Clim <- c(0,4);
    cex <- w / mean(w, na.rm=TRUE) + 1/2;
    plot(C1C2, cex=cex, xlim=Clim, ylim=Clim);
    lines(C1C2);
    abline(a=0, b=1, lty=3);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Transform (C1,C2) -> (dC1, dC2) -> (a,b) line space
  #
  # Lines can be represented as:
  #  (1) Two points (x0,y0) and (x1,y1).
  #  (2) Intercept and slope (a,b).
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Represent change points in line space");
  # Change-point weights
  cpw <- w[1:(length(w)-1)] + w[2:length(w)];

  # Alternative (1a): y = a + b*x
  # Slope    : b = dy/dx, where dy=(y1-y0) and dx=(x1-x0)
  # Intercept: a = y - b*x = y0 - b*x0

  # Alternative (1b): y = a + b*x
  # Slope    : beta = arctan(dx/dy)
  #              s.t. beta=0 <=> vertical, beta=pi/2 <=> horizontal
  # Intercept: a = y - b*x = y0 - b*x0

  AB <- C1C2toAB(X);

  verbose && cat(verbose, "(C1,C2) -> (A,B):");
  verbose && print(verbose, AB);

  verbose && cat(verbose, "Slopes in (C1,C2) for all change points:");
  b <- AB[,"b"];
  verbose && print(verbose, b);

  # Keep only finite data points
  ok <- (is.finite(b) & is.finite(cpw));
  bT <- b[ok];
  verbose && cat(verbose, "Finite slopes in (C1,C2):");
  cpwT <- cpw[ok];
  cpwT <- cpwT / sum(cpwT);
  verbose && print(verbose, cbind(bT=bT,wT=cpwT));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Adjust modes in {b} to their expect locations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Estimate density
  jitter <- 0;
#  jitter <- rnorm(length(bT), mean=0, sd=0.1);
#  bT <- bT + jitter;
  adjust <- 0.8;

  if (debugPlot) {
    bLim <- c(-pi/2, pi/2);
    xlim <- bLim + 0.1*bLim;
    # Draw density before adjustment
    plotDensity(bT, weights=cpwT, lwd=2, from=bLim[1], to=bLim[2], xlim=xlim, ylim=c(0,2), xlab="slopes", adjust=adjust);
#    plotDensity(bT, adjust=2);
    abline(v=0, lty=3); # Vertical lines
    abline(v=c(-pi/2,pi/2), lty=1); # Horizontal lines
    abline(v=c(-pi/4,pi/4), lty=3); # Diagonal lines
    x0 <- par("usr")[1]; x1 <- par("usr")[2];
    y0 <- par("usr")[3]; y1 <- par("usr")[4];
    dx <- (x1-x0); dy <- (y1-y0);
    # Annotate slopes
    yy <- c(0.85,0.90)*dy;
    xx <- 0 + c(0,0);
    arrows(x0=xx[1], x1=xx[2], y0=yy[1], y1=yy[2], code=2, lwd=2, length=0.05, col="blue");
    xx <- -pi/2 + c(0,-0.05*dx) + 0.05/2*dx;
    arrows(x0=xx[1], x1=xx[2], y0=yy[1], y1=yy[1], code=2, lwd=2, length=0.05, col="blue");
    xx <- +pi/2 + c(0,+0.05*dx) - 0.05/2*dx;
    arrows(x0=xx[1], x1=xx[2], y0=yy[1], y1=yy[1], code=2, lwd=2, length=0.05, col="blue");
    xx <- -pi/4 + c(0,-0.05*dx) + 0.05/2*dx;
    arrows(x0=xx[1], x1=xx[2], y0=yy[1], y1=yy[2], code=2, lwd=2, length=0.05, col="blue");
    xx <- +pi/4 + c(0,+0.05*dx) - 0.05/2*dx;
    arrows(x0=xx[1], x1=xx[2], y0=yy[1], y1=yy[2], code=2, lwd=2, length=0.05, col="blue");
  }

  # Find modes in {b}
  fp <- findPeaksAndValleys(bT, weights=cpwT, from=bLim[1], to=bLim[2], tol=0.05, adjust=adjust);
  verbose && cat(verbose, "Peaks and valleys:");
  verbose && print(verbose, fp);
  type <- NULL; rm(type); # To please R CMD check
  pfp <- subset(fp, type == "peak");
  xPeaks <- pfp[,"x"];
  verbose && cat(verbose, "Slopes for all modes: ", 
                        paste(sprintf("%.3f", xPeaks), collapse=", "));

  # Call modes
  # Alt 1
  expected <- c(-1/2,-1/4,0,+1/4,+1/2)*pi;
  xs <- xPeaks;
  calls <- sapply(xs, FUN=function(x) {
    dist <- abs(x - expected);
    which.min(dist);
  });
  names(xs) <- expected[calls];
  # Alt 2:
  # It is probably better to call the strongest peaks first for which
  # we have more confidence, and then call the other relative to those.
  # /HB 2010-09-19
  expected <- c(-1/2,-1/4,0,+1/4,+1/2)*pi;
  xd <- pfp[,c("x", "density"),drop=FALSE];
#  xd[,"x"] <- xd[,"x"];
  xd <- xd[order(xd[,"density"], decreasing=TRUE),,drop=FALSE];
  calls <- rep(NA, times=nrow(xd));
  expectedLeft <- expected;
  for (kk in seq(length=nrow(xd))) {
    # Mode #kk
    x <- xd[kk,"x"];
    dx <- abs(x - expectedLeft);
    call <- which.min(dx);
    expectedLeft[call] <- NA;
    calls[kk] <- call;
  } # for (kk ...)
  xde <- cbind(xd, expected=expected[calls]);
  if (debugPlot) {
    # Highlight and annotate modes
    d <- apply(xde[,c("x", "density")], MARGIN=1, FUN=function(xd) {
      lines(x=rep(xd[1], times=2), y=c(y0,xd[2]), lwd=2);
    });
    text(xde[,c("x","density")], labels=seq(length=nrow(xde)), adj=c(0.5,-1));

    # Annotate deviance from estimated slopes
    yT <- y0 + 0.05*dy;
    apply(xde[,c("x", "expected")], MARGIN=1, FUN=function(xy) {
      arrows(x0=xy[1], x1=xy[2], y0=yT, y1=yT, code=2, lwd=2, length=0.08, col="red");
    });
  }

  # Fit backtransform
  X <- as.matrix(xde[,c("x","expected")]);
  w <- xde[,"density",drop=TRUE];
  w <- w / sum(w);
  n <- nrow(X);
  if (n == 1) {
    # Global shift
    shift <- X[,"expected"] - X[,"x"];
    fitB <- list(predictY = function(x) x+shift);
  } else if (n == 2) {
    x <- X[,1];
    y <- X[,2];
    fitB <- lm(y~x, weights=w);
    fitB$predictY <- function(x) {
      predict(fitB, newdata=data.frame(x=x));
    }
  } else {
  #  fitB <- fitXYCurve(X, weights=w, method="loess"); # ?!?
    fitB <- fitXYCurve(X, method="lowess");
  }

  # Backtransform
  bTHat <- fitB$predictY(bT);

  if (debugPlot) {
    ok <- (is.finite(bTHat) & is.finite(cpwT));
    # Draw density after adjustment
    plotDensity(bTHat[ok], weights=cpwT[ok], from=bLim[1], to=bLim[2], adjust=adjust, col="red", lty=3, add=TRUE);
  }

  if (debugPlot) {
    plot(X, cex=2*w+1/2, xlim=bLim, ylim=bLim, xlab="slopes", ylab="expected");
    abline(a=0, b=1, lty=3, col="gray");
    x <- seq(from=bLim[1],to=bLim[2], by=0.01);
    y <- fitB$predictY(x);
    lines(x,y, lwd=2);
  }

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Orthogonalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Orthogonalize the (C1,C2) space by adjusting slopes to their expected locations");
  # Backtransform
  bHat <- fitB$predictY(b);
  verbose && print(verbose, b);
  verbose && print(verbose, bHat);

  # Verticalization
  #  b = atan(dC1C2[,1]/dC1C2[,2]);
  #  tan(b) = dC1/dC2;
  #  dC2*tan(b) = dC1;
  abHat <- AB;
  abHat[,2] <- bHat;
  C1C2o <- ABtoC1C2(abHat, C1C2=C1C2);
  verbose && print(verbose, C1C2o-C1C2);

  if (debugPlot) {
    Clim <- c(0,4);
    plot(C1C2o, cex=cex, xlim=Clim, ylim=Clim);
    lines(C1C2o);
    abline(a=0, b=1, lty=3);
  }

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

  verbose && exit(verbose);

  fitO;
}) # orthogonalizeC1C2()


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
setMethodS3("deShearC1C2", "PairedPSCBS", function(fit, weightFlavor=c("min", "sum"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'weightFlavor':
  weightFlavor <- match.arg(weightFlavor);
  
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
  X <- extractMinorMajorCNs(fit);

  # Number of TCN and DH data points
  counts <- X[,3:4, drop=FALSE];

  # Region weights from DH counts
  w <- counts[,2];
  w <- sqrt(w);
  w <- w / sum(w, na.rm=TRUE);

  ## Change-point weights
  if (weightFlavor=="min") {
    ## smallest of the two flanking (DH) counts 
    cpw <- cbind(w[1:(length(w)-1)], w[2:length(w)]);
    cpw <- rowMins(cpw);
    cpw <- sqrt(cpw);
  } else if (weightFlavor=="sum") {
    ## sum of region weights
    cpw <- w[1:(length(w)-1)] + w[2:length(w)];
  }
  
  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];
  D <- colDiffs(C1C2);
  alpha <- atan(D[, 2]/D[, 1]);

  verbose && enter(verbose, "Represent change points as angles");
  verbose && str(verbose, alpha);

  # Keep only finite data points
  ok <- (is.finite(alpha) & is.finite(cpw));
  alphaT <- alpha[ok];
  verbose && cat(verbose, "Finite slopes in (C1,C2):");
  cpwT <- cpw[ok];
  cpwT <- cpwT / sum(cpwT);
  verbose && print(verbose, cbind(alphaT=alphaT,wT=cpwT));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Adjust modes in {alpha} to their expect locations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ## The signal is pi-periodic.
  ## we are looking at from -pi/2 to pi/2.
  ## we expect a peak near -pi/2 (or pi/2...)
  ## in order to estimate it correctly, transform the signal so that it is
  ## in -pi/2-pi/8, pi/2-pi/8
  aa <- alphaT;
  ww <- which(alphaT > (pi/2+pi/4)/2) ## half way to both expected peaks
  aa[ww] <- aa[ww]-pi;
  rg <- range(aa, na.rm=TRUE);

  fp <- findPeaksAndValleys(aa, weights=cpwT, from=rg[1], to=rg[2], ...);
  verbose && cat(verbose, "Peaks and valleys:");
  verbose && print(verbose, fp);
  type <- NULL; rm(type); # To please R CMD check
  pfp <- subset(fp, type == "peak");

  xPeaks <- pfp[,"x"];
  
  ## Call modes
  expected <- c(-1/2,-1/4,0,+1/4,+1/2)*pi;
  xs <- xPeaks;
  calls <- sapply(xs, FUN=function(x) {
    dist <- abs(x - expected);
    which.min(dist);
  });
  names(xs) <- expected[calls];
  ## TODO: more robust ?

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Orthogonalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Remove shearing from the (C1,C2) space by independent linear adjustment in x and y");
  ## Fit backtransform
  H <- function(xy, thetaXY) {
    sx <- tan(thetaXY[1]+pi/2);
    sy <- tan(thetaXY[2]+pi/2);
    cbind(xy[, 1]+xy[, 2]*sx, xy[, 1]*sy + xy[, 2]);
  }

  ## Backtransform
  ## Use -pi/2 and 0 to correct for shearing
  wx <- calls[match(-pi/2, expected)];
  wy <- calls[match(0, expected)];
  
  tX <- xPeaks[wx]
  tY <- pi/2 - xPeaks[wy]

  C1C2o <- H(C1C2, c(tX, tY));
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

  verbose && exit(verbose);

  fitO;
}) # deShearC1C2()

##############################################################################
# HISTORY
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

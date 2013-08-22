# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For (C1,C2) data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("estimateChangePointFlips", "matrix", function(C, ...) {
  estimateChangePointFlips(as.data.frame(C), ...);
}, protected=TRUE)


# TODO: Return/incorporate information on flanking segment lengths.
#       This is part of 'data', but not used and not returned right
#       now.  Decision is to ignore that now.  Maybe we'll later
#       use bootstrapping or similar to estimate confidence scores
#       which in turn may affect the call. /HB+PN 2013-08-15 in Banff.
setMethodS3("estimateChangePointFlips", "data.frame", function(data, ...) {
  n <- nrow(data);
  nbrOfCPs <- n - 1L;
  alpha <- matrix(NA_real_, nrow=nbrOfCPs, ncol=2L);
  colnames(alpha) <- c("alpha0", "alpha1");

  # (C1,C2)
  C <- as.matrix(data[,1:2,drop=FALSE]);

  # For each change point
  for (ii in seq_len(nbrOfCPs)) {
    rows <- ii + 0:1;

    # alpha0: angle without flip
    C0 <- C[rows,];
    D0 <- C0[2,] - C0[1,];
    alpha0 <- atan(D0[1]/D0[2]);

    # alpha1: angle with flip
    C1 <- C0; C1[2,] <- C1[2,2:1];
    D1 <- C1[2,] - C1[1,];
    alpha1 <- atan(D1[1]/D1[2]);

    alpha[ii,1L] <- alpha0;
    alpha[ii,2L] <- alpha1;
  } # for (ii ...)

  betaA <- abs(alpha - pi/4);
  colnames(betaA) <- c("betaA0", "betaA1");
  betaB <- abs(alpha + pi/4);
  colnames(betaB) <- c("betaB0", "betaB1");
  gammaA <- pmin(betaA[,1], betaB[,1]);
  gammaB <- pmin(betaA[,2], betaB[,2]);
  call <- (gammaB > gammaA);
  cbind(alpha, betaA=betaA, betaB=betaB, gammaA=gammaA, gammaB=gammaB, call=call);
}, protected=TRUE)


setMethodS3("callChangePointFlips", "matrix", function(C, ...) {
  callChangePointFlips(as.data.frame(C, ...));
}, protected=TRUE)

setMethodS3("callChangePointFlips", "data.frame", function(C, ...) {
  theta <- estimateChangePointFlips(C, ...);
  flip <- as.logical(theta[,"call"]);
  which(flip) + 1L;
}, protected=TRUE)

setMethodS3("callChangePointFlips", "PSCBS", function(fit, ...) {
  C <- extractC1C2(fit);
  flips <- callChangePointFlips(C);
  flips;
})


setMethodS3("adjustAB", "PSCBS", function(fit, ...) {
  segs <- getSegments(fit);
  segs$dhMean[segs$abCall] <- 0;
  fit$output <- segs;
  params <- fit$params;
  params$abAdjusted <- TRUE;
  fit$params <- params;
  fit;
}, protected=TRUE)


setMethodS3("flipChangePoints", "PSCBS", function(fit, flips=callChangePointFlips(fit), ...) {
  # Argument 'flips':
  flips <- Arguments$getIndices(flips, max=nbrOfSegments(fit, splitters=TRUE));

  # Update segment data
  if (length(flips) > 0L) {
    segs <- getSegments(fit);
    fields <- colnames(segs);

    # (a) Update (C1,C2) fields - swap (C1,C2)
    c1Names <- grep("^c1", fields, value=TRUE);
    c2Names <- gsub("^c1", "c2", c1Names);
    stopifnot(length(c2Names) == length(c1Names));

    for (kk in seq_along(c1Names)) {
      cols <- c(c1Names[kk], c2Names[kk]);
      C <- segs[,cols];
      C[flips,] <- C[flips,2:1];
      segs[,cols] <- C;
    }

##    # (b) Update DH fields - swap signs
##    dhNames <- grep("^dh", fields, value=TRUE);
##    dhNames <- setdiff(dhNames, c("dhId", "dhStart", "dhEnd", "dhNbrOfLoci"));
##    for (name in dhNames) {
##      segs[[name]] <- -segs[[name]];
##    }

    segs$c1c2Swap <- FALSE;
    segs$c1c2Swap[flips] <- TRUE;

    fit$output <- segs;
  }

  # Update parameter estimates
  params <- fit$params;
  params$flipChangePoints <- list(
    flips = flips
  );
  fit$params <- params;

  fit;
}, protected=TRUE)


setMethodS3("deshearC1C2", "PSCBS", function(fit, S=0, ...) {
  # Argument 'S':
  S <- Arguments$getNumeric(S);

  # Deshear (C1,C2) which is dual to adjusting DH.
  segs <- getSegments(fit, splitters=TRUE);
  fields <- colnames(segs);

  # Update DH fields
  dhNames <- grep("^dh", fields, value=TRUE);
  dhNames <- setdiff(dhNames, c("dhId", "dhStart", "dhEnd", "dhNbrOfLoci"));
  for (name in dhNames) {
    rho <- segs[[name]];;
    rho <- (1+S)/(1-S)*rho;
    segs[[name]] <- rho;
  }

  fit$output <- segs;

  # Update parameters
  params <- fit$params;
  params$deshearC1C2 <- list(S=S);
  fit$params <- params;

  fit;
}, protected=TRUE) # deshearC1C2()


setMethodS3("estimateDeshearingParameter_20130815", "PSCBS", function(fit, ...) {
  D <- extractDeltaC1C2(fit);

  ## == atan(dY/dX)
  alpha <- atan(D[,2L]/D[,1L]);

  ## this is tan(alpha+pi/4)
  x <- (D[,2L]-D[,1L])/(D[,2L]+D[,1L]);

  ## absolute value should be safe because we are looking for a mode
  ## close to pi/4.
  beta <- abs(atan(x)+pi/4);

  fitD <- findPeaksAndValleys(beta);

  type <- NULL; rm(list="type");  # To please R CMD check
  peaks <- subset(fitD, type=="peak");

  ## the closest peak to 0
  ww <- which.min(abs(peaks$x));
  xx <- peaks$x[ww];

  ## the deshearing parameter
  tan(xx);
}, protected=TRUE) # estimateDeshearingParameter_20130815()

setMethodS3("estimateDeshearingParameter", "PSCBS", function(fit, ...) {
  D <- extractDeltaC1C2(fit);
  alpha <- atan(D[,2L]/D[,1L]);
  alpha <- 360/(2*pi)*alpha;          ## converting to degrees
  beta <- (alpha + 45) %% 180 - 45;   ## mapping -pi/2 mode to pi/2
  beta <- 45 - abs(beta-45);          ## shearing is supposed to be
                                      ## symmetric wrt 45

  beta <- na.omit(beta);
  fitD <- findPeaksAndValleys(beta);
  type <- NULL; rm(list = "type");    ## To please R CMD check
  peaks <- subset(fitD, type == "peak");
  ww <- which.min(abs(peaks$x));      ## closest peak to 0
  xx <- peaks$x[ww];
  atan(xx*(2*pi)/360);
}, protected=TRUE) # estimateDeshearingParameter()

##############################################################################
# HISTORY
# 2013-08-21
# o Updated estimateDeshearingParameter().
# o Added PN's sketch on estimateDeshearingParameter().
# 2013-08-15
# o Added to aroma.cn.
# 2013-08-09
# o Created.
##############################################################################

setMethodS3("extractPolarDeltaC1C2", "PairedPSCBS", function(..., modulo=c("2pi", "pi")) {
  # Argument 'modulo':
  modulo <- match.arg(modulo);

  data <- extractDeltaC1C2(...);
  radius <- sqrt(data[,2]^2 + data[,1]^2);
  alpha <- atan(data[,2]/data[,1]);

  # Make angular into [0,2*pi]
  alpha <- alpha + base::pi/2;

  if (modulo == "pi") {
    alpha <- alpha %% base::pi;
  }

  data[,1] <- radius;
  data[,2] <- alpha;

  data;
})


setMethodS3("extractDeltaC1C2", "PairedPSCBS", function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  X <- extractC1C2(...);

  # (C1,C2)
  C1C2 <- X[,1:2,drop=FALSE];

  # Number of TCN and DH data points
  counts <- X[,3:4, drop=FALSE];

  # Not needed anymore
  rm(X);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate (dC1,dC2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (dC1, dC2)
  dC1C2 <- matrixStats::colDiffs(C1C2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Change-point weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Region weights from DH counts
  w <- counts[,2,drop=TRUE];
  w <- sqrt(w);
  w <- w / sum(w, na.rm=TRUE);

  # (a) Smallest of the two flanking (DH) counts 
  cpw <- cbind(w[1:(length(w)-1)], w[2:length(w)]);
  cpw <- rowMins(cpw, na.rm=TRUE);
  cpw[is.infinite(cpw)] <- NA;
  cpw <- sqrt(cpw);
  cpwMin <- cpw / sum(cpw, na.rm=TRUE);

  # (b) Sum of region weights
  cpw <- w[1:(length(w)-1)] + w[2:length(w)];
  cpwAvg <- cpw / sum(cpw, na.rm=TRUE);

  cbind(dC1=dC1C2[,1], dC2=dC1C2[,2], wMin=cpwMin, wAvg=cpwAvg);
})



setMethodS3("extractC1C2", "list", function(fitList, ...) {
  c1c2List <- lapply(fitList, FUN=function(fit) {
    extractC1C2(fit);
  });
  
  # Append NAs between chromosomes
  c1c2TList <- lapply(c1c2List, FUN=function(c1c2) {
    rbind(c1c2, NA);
  })
  
  c1c2 <- Reduce(rbind, c1c2TList);

  w <- sqrt(c1c2[,4]);
  w <- w / sum(w, na.rm=TRUE);
  w <- w / mean(w, na.rm=TRUE);
  c1c2 <- cbind(c1c2, w=w);

  class(c1c2) <- unique(c("C1C2", class(c1c2)));

  c1c2;
});

setMethodS3("plot", "C1C2", function(x, xlim=c(0,4), ylim=xlim, xlab=expression(C[1]), ylab=expression(C[2]), ...) {
  plot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
  abline(a=0, b=1, lty=3);
  grid();
  points(x, ...);
})

setMethodS3("points", "C1C2", function(x, cex=sqrt(x[,"w"])+1/8, ...) {
  NextMethod("points", x[,c("C1","C2"),drop=FALSE], cex=cex, ...);
})


setMethodS3("fitLoess2D", "C1C2", function(X, Y, ...) {
  fit <- fitLoessKD(X=X[,1:2,drop=FALSE], Y=Y[,1:2,drop=FALSE]);
  class(fit) <- c("Loess2DFit", class(fit));
  fit;
})

setMethodS3("normalizeLoess2D", "C1C2", function(X, ...) {
  XN <- X;
  XN[,1:2] <- normalizeLoessKD(X[,1:2], ...);
  XN;
}) # normalizeLoess2D()


##############################################################################
# HISTORY
# 2010-10-05 [HB]
# o Added extractPolarDeltaC1C2(), which returns [0,pi] angles and lengths.
# o Now extractDeltaC1C2() returns two types of change-point weights.
# 2010-09-19 [HB]
# o Created.
##############################################################################

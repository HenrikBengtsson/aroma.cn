testEqualityTcnByT <- function(dataL, dataR, alpha=0.02, ...) {
  nL <- length(dataL);
  nR <- length(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }
  print(c(nL=nL, nR=nR));

  # Test:
  #  H0: muL == muR
  #  H1: muL != muR
  fit <- t.test(dataL, dataR, paired=FALSE, var.equal=TRUE, 
                                            alternative="two.sided");
  t <- fit$statistic;
  p <- fit$p.value;
  fit$isSignificant <- (p < alpha);
  isEqual <- (!fit$isSignificant);
  attr(isEqual, "fit") <- fit;

  isEqual;
} # testEqualityTcnByT()


testEqualityTcnByMean <- function(dataL, dataR, alpha=0.02, ...) {
  nL <- length(dataL);
  nR <- length(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }

  muL <- mean(dataL, na.rm=TRUE);
  muR <- mean(dataR, na.rm=TRUE);
  delta <- abs(muL - muR);

  isEqual <- (delta <= alpha);

  isEqual;
} # testEqualityTcnByMean()


testEqualityC1C2ByMean <- function(dataL, dataR, alpha=0.02, ...) {
  # Sanity checks
  stopifnot(is.matrix(dataL));
  stopifnot(is.matrix(dataR));

  nL <- nrow(dataL);
  nR <- nrow(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }

  muL <- colMeans(dataL, na.rm=TRUE);
  muR <- colMeans(dataR, na.rm=TRUE);
  delta <- abs(muL - muR);

  # Sanity checks
  stopifnot(length(delta) == 2);

  # Eucledian distance
  delta <- sqrt(sum(delta^2, na.rm=TRUE));

  isEqual <- (delta <= alpha);

  isEqual;
} # testEqualityC1C2ByMean()


testEqualityC1C2ByChiSquare <- function(dataL, dataR, alpha=0.02, ...) {
  # Sanity checks
  stopifnot(is.matrix(dataL));
  stopifnot(is.matrix(dataR));

  nL <- nrow(dataL);
  nR <- nrow(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }

  stop("Not yet implemented!");

  # Test:
  #  H0: muL == muR
  #  H1: muL != muR
##  fit <- t.test(dataL, dataR, paired=FALSE, var.equal=TRUE, 
##                                            alternative="two.sided");

  t <- fit$statistic;
  p <- fit$p.value;
  fit$isSignificant <- (p < alpha);
  isEqual <- (!fit$isSignificant);
  attr(isEqual, "fit") <- fit;

  isEqual;
} # testEqualityC1C2ByChiSquare()



############################################################################
# HISTORY:
# 2011-01-18
# o Added testEqualityC1C2ByMean().
# 2011-01-12
# o Added testEqualityTCN().
# o Created.
############################################################################

setMethodS3("findPeaksAndValleys", "numeric", function(x, ..., tol=0.01) {
  # Argument 'tol':
  tol <- Arguments$getDouble(tol, range=c(.Machine$double.eps, Inf));

  d <- density(x, ...);
  delta <- diff(d$y);
  n <- length(delta);

  isPeak <- (delta[-n] > 0 & delta[-1] < 0);
  isValley <- (delta[-n] < 0 & delta[-1] > 0);
  isPeakOrValley <- (isPeak | isValley);

  idxs <- which(isPeakOrValley);
  types <- c("valley", "peak")[isPeak[idxs]+1];
  names(idxs) <- types;

  x <- d$x[idxs];
  y <- d$y[idxs];
  res <- data.frame(type=types, x=x, density=y);
  res <- subset(res, density >= tol);
  res;
}) # findPeaksAndValleys()

############################################################################
# HISTORY:
# 2009-03-06 [HB]
# o Created for doing quick naive genotyping of some TCGA normal samples in
#   order to highlight the centers of the clouds in a tumor-normal fracB
#   scatter plots.
############################################################################

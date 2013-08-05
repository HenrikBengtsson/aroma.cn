setMethodS3("makeSmoothSplinePredict", "numeric", function(x, y, df=5, ...) {
  .Defunct("makeSmoothSplinePredict", msg="Internal function makeSmoothSplinePredict() is deprecated, because it is no longer used.  It will eventually be removed from the package.");

  # Argument 'y':
  if (length(x) != length(y)) {
    throw("Argument 'y' is of a different length than 'x'");
  }

  # Identify finite (x,y) pairs
  ok <- which(is.finite(x) & is.finite(y));
  x <- x[ok];
  y <- y[ok];
  # Not needed anymore
  ok <- NULL;

  specs <- list(
    xRange = range(x),
    yRange = range(y)
  );

  # Fit smooth function
  fit <- smooth.spline(x,y, df=df, ...);
  # Not needed anymore
  x <- y <- NULL;

  # Create predict() function that handles missing values.
  predFcn <- function(x, ...) {
    yPred <- rep(as.double(NA), times=length(x));
    ok <- which(!is.na(x));  # Allows for -/+Inf:s though
    yPred[ok] <- predict(fit, x=x[ok], ...)$y;
    yPred;
  } # predFcn()

  attr(predFcn, "fit") <- fit;
  attr(predFcn, "specs") <- specs;

  predFcn;
}, private=TRUE, deprecated=TRUE) # makeSmoothSplinePredict()



###########################################################################
# HISTORY:
# 2013-08-04
# o CLEAN UP: Defuncted makeSmoothSplinePredict().
# 2009-02-08
# o CLEAN UP: Deprecated makeSmoothSplinePredict(), because it is not used.
# 2008-10-08
# o Removed implementation for data.frame:s.
# 2008-05-27
# o Created.  Will probably end up in aroma.light.
###########################################################################

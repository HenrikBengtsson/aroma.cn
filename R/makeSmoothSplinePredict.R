setMethodS3("makeSmoothSplinePredict", "numeric", function(x, y, df=5, ...) {
  # Argument 'y':
  if (length(x) != length(y)) {
    throw("Argument 'y' is of a different length than 'x'");
  }

  # Identify finite (x,y) pairs
  ok <- which(is.finite(x) & is.finite(y));
  x <- x[ok];
  y <- y[ok];
  rm(ok);

  specs <- list(
    xRange = range(x),
    yRange = range(y)
  );

  # Fit smooth function
  fit <- smooth.spline(x,y, df=df, ...);
  rm(x,y);

  # Create predict() function that handles missing values.
  predFcn <- function(x, ...) {
    yPred <- rep(as.double(NA), length(x));
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
# 2009-02-08
# o CLEAN UP: Deprecated makeSmoothSplinePredict(), because it is not used.
# 2008-10-08
# o Removed implementation for data.frame:s.
# 2008-05-27
# o Created.  Will probably end up in aroma.light.
###########################################################################

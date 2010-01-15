setMethodS3("getBacktransforms", "PrincipalCurve", function(fit, dimensions=NULL, targetDimension=1, range=NULL, length.out=100, ...) {
  # Argument 'dimensions':
  s <- fit$s;
  ndim <- ncol(s);
  if (is.null(dimensions)) { 
    dimensions <- seq(length=ndim);
  }
  dimensions <- Arguments$getIndices(dimensions, max=ndim);

  # Argument 'range':
  if (is.null(range)) {
    range <- range(s, na.rm=TRUE);  
  }
  range <- Arguments$getDoubles(range, length=c(2,2), disallow=c("NA", "NaN"));


  y <- seq(from=range[1], to=range[2], length.out=length.out);

  naValue <- as.double(NA);
  dim <- c(length(y), 2, length(dimensions));
  XY <- array(naValue, dim=dim);

  for (kk in seq(along=dimensions)) {
    dim <- dimensions[kk];
    yN <- backtransformPrincipalCurve(y, fit=fit, dimensions=dim,
                                     targetDimension=targetDimension);
    yN <- yN[,1,drop=TRUE];
    xy <- cbind(y, yN);

    XY[,,kk] <- xy;
  } # for (kk ...)

  XY;
})

setMethodS3("plotBacktransforms", "PrincipalCurve", function(fit, ..., xlim=c(-3,3), ylim=xlim, xlab="y", ylab="y*") {
  XY <- getBacktransforms(fit, ...);

  if (is.null(xlim)) {
    range(XY, na.rm=FALSE);
  }
  if (is.null(ylim)) {
    ylim <- xlim;
  }

  ndim <- dim(XY)[3];
  subplots(ndim);
  for (kk in seq(length=ndim)) {
    xy <- XY[,,kk];
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
    abline(a=0, b=1, lty=3);
    lines(xy, col="red", lwd=2);
  }

  invisible(XY);
})


############################################################################
# HISTORY:
# 2010-01-14
# o Added getBacktransforms() and plotBacktransforms().
# o Created.
############################################################################

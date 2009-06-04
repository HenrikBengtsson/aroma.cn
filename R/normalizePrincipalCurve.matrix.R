setMethodS3("normalizePrincipalCurve", "matrix", function(x, ..., center=TRUE, returnFit=FALSE) {
  # Fit principal curve
  fit <- fitPrincipalCurve(x, ...);

  # Flip direction of 'lambda'?
  rho <- cor(fit$lambda, x[,1], use="complete.obs");
  flip <- (rho < 0);

  dx <- (fit$s - x);
  xN <- fit$lambda + dx;
  if (flip) {
    xN <- -xN;
  }

  if (center) {
    # Same center for each column
    for (cc in seq(length=ncol(x))) {
      mu <- median(x[,cc], na.rm=TRUE);
      muN <- median(xN[,cc], na.rm=TRUE);
      xN[,cc] <- xN[,cc] - (muN-mu);
    }
  }

  # Return fit?
  if (returnFit)
    attr(xN, "fit") <- fit;

  xN;
}) # normalizePrincipalCurve()



###########################################################################
# HISTORY:
# 2008-10-08
# o Removed implementation for data.frame:s.
# 2008-05-27
# o Added normalizePrincipalCurve().
# o Created.  Will probably end up in aroma.light.
###########################################################################

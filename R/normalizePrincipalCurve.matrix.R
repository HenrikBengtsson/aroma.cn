#########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizePrincipalCurve
#
# \encoding{latin1}
#
# @title "Fit a principal curve in K dimensions"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{x}{An NxK @matrix (K>=2) where the columns represent the dimension.}
#  \item{...}{Additional arguments passed to
#      @see "aroma.light::fitPrincipalCurve" used for fitting the model.}
#  \item{center}{If @TRUE, normalized data is centered such that the median
#      signal in each dimension is at zero.}
#  \item{returnFit}{If @TRUE, the fitted principal curve parameters are
#      returned as an attribute.}
# }
#
# \value{
#   Returns an NxK @matrix.
# }
#
# @author
#
# \references{
#   [1] Hastie, T. and Stuetzle, W, \emph{Principal Curves}, JASA, 1989.
# }
#
# \seealso{
#   @see "aroma.light::fitPrincipalCurve" and
#   @see "aroma.light::backtransformPrincipalCurve".
# }
#*/#########################################################################  
setMethodS3("normalizePrincipalCurve", "matrix", function(x, ..., center=TRUE, returnFit=FALSE) {
  # fitPrincipalCurve()
  require("aroma.light") || throw("Package not loaded: aroma.light"); 

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
# 2012-04-16
# o Now normalizePrincipalCurve() explicitly require aroma.light.
# 2012-03-30
# o Added Rdoc comments.
# 2008-10-08
# o Removed implementation for data.frame:s.
# 2008-05-27
# o Added normalizePrincipalCurve().
# o Created.  Will probably end up in aroma.light.
###########################################################################

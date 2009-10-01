###########################################################################/**
# @set "class=list"
# @RdocMethod normalizeDifferencesToAverage
# @alias normalizeDifferencesToAverage
#
# @title "Rescales channel vectors to get the same average"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{A @numeric @list of length K.}
#   \item{baseline}{An @integer in [1,K] specifying which channel should be
#      the baseline.  The baseline channel will be almost unchanged.
#      If @NULL, the channels will be shifted towards median of them all.}
#   \item{FUN}{A @function for calculating the average of one channel.}
#   \item{...}{Additional arguments passed to the \code{avg} @function.}
# }
#
# \value{
#  Returns a normalized @list of length K.
# }
#
# @examples "../incl/normalizeDifferencesToAverage.Rex"
#
# @author
#*/########################################################################### 
setMethodS3("normalizeDifferencesToAverage", "list", function(x, baseline=1, FUN=median, ...) {
  # Argument 'x':
  if (!is.list(x)) {
    throw("Argument 'x' is not a list: ", class(x)[1]);
  }

  # Number dimensions
  ndim <- length(x);

  # Argument 'baseline':
  if (!is.null(baseline)) {
    baseline <- Arguments$getIndex(baseline, range=c(1,ndim));
  }


  # Calculate the overall average level for each dimension
  mus <- sapply(x, FUN=function(y) {
    y <- y[is.finite(y)];
    FUN(y);
  });

  # Estimate the overall target level
  if (is.null(baseline)) {
    targetMu <- mus[baseline];
  } else {
    targetMu <- median(mus, na.rm=TRUE);
  }

  # Calculate the required overall shifts per dimension
  deltas <- mus - targetMu;

  # Shift all dimensions so that all have the same overall average
  xN <- mapply(x, as.list(deltas), FUN=function(y, dy) {
    y <- y - dy;
    list(y);
  });

  # Return estimated parameters
  fit <- list(mus=mus, baseline=baseline, targetMu=targetMu, deltas=deltas);
  attr(xN, "fit") <- fit;

  xN;
})


############################################################################
# HISTORY:
# 2009-09-30
# o Created from the source of an aroma.tcga vignette from May 2009.
############################################################################

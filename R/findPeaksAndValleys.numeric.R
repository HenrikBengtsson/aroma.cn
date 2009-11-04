###########################################################################/**
# @set "class=numeric"
# @RdocMethod findPeaksAndValleys
# @alias findPeaksAndValleys
#
# @title "Finds extreme points in the empirical density estimated from data"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{x}{A @numeric @vector containing data points.}
#  \item{na.rm}{If @TRUE, missing values are dropped, otherwise not.}
#  \item{...}{Arguments passed to @see "stats::density".}
#  \item{tol}{A non-negative @numeric threshold specifying the minimum
#    density at the extreme point in order to accept it.}
# }
#
# \value{
#   Returns a @dataframe of "peaks" and "valleys" filtered by \code{tol}.
# }
#
# @examples "..\incl\findPeaksAndValleys.Rex"
#
# @author
#*/########################################################################### 
setMethodS3("findPeaksAndValleys", "numeric", function(x, na.rm=TRUE, ..., tol=0) {
  # Argument 'na.rm':
  na.rm <- Arguments$getLogical(na.rm);

  # Argument 'tol':
  tol <- Arguments$getDouble(tol, range=c(0, Inf));

  d <- density(x, na.rm=na.rm, ...);
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

  # Filter valleys by density?
  if (tol > 0) {
    res <- subset(res, density >= tol);
  }

  res;
}) # findPeaksAndValleys()

############################################################################
# HISTORY:
# 2009-11-03 [HB]
# o Added Rdoc comments with an example().
# 2009-03-06 [HB]
# o Created for doing quick naive genotyping of some TCGA normal samples in
#   order to highlight the centers of the clouds in a tumor-normal fracB
#   scatter plots.
############################################################################

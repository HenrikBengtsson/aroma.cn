###########################################################################/**
# @set "class=numeric"
# @RdocMethod callNaiveGenotypes
# @alias callNaiveGenotypes
#
# @title "Calls genotypes in a normal sample"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{y}{A @numeric @vector of length J containing allele B fractions
#    for a normal sample.}
#  \item{cn}{An optional @numeric @vector of length J specifying the true
#    total copy number in \eqn{\{0,1,2,NA\}} at each locus.  This can be 
#    used to specify which loci are diploid and which are not, e.g. 
#    autosomal and sex chromosome copy numbers.}
#  \item{flavor}{A @character string specifying the type of algorithm used.}
#  \item{adjust}{A postive @double specifying the amount smoothing for
#    the empirical density estimator.}
#  \item{...}{Additional arguments passed to @see "findPeaksAndValleys".}
#  \item{censorAt}{A @double @vector of length two specifying the range
#    for which values are considered finite.  Values below (above) this 
#    range are treated as -@Inf (+@inf).}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a @numeric @vector of length J containing the genotype calls
#   in allele B fraction space, that is, in [0,1] where 1/2 corresponds
#   to a heterozygous call, and 0 and 1 corresponds to homozygous A 
#   and B, respectively.
#   Non called genotypes have value @NA.
# }
#
# @examples "..\incl\callNaiveGenotypes.Rex"
#
# \section{Missing and non-finite values}{
#   A missing value always gives a missing (@NA) genotype call.
#   Negative infinity (-@Inf) always gives genotype call 0.
#   Positive infinity (+@Inf) alwsys gives genotype call 1.
# }
#
# @author
#
# \seealso{
#   Internally @see "findPeaksAndValleys" is used to identify the thresholds.
# }
#*/########################################################################### 
setMethodS3("callNaiveGenotypes", "numeric", function(y, cn=rep(2L, length(y)), flavor=c("density"), adjust=1.5, ..., censorAt=c(-0.5,+1.5), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  J <- length(y);
  y <- as.double(y);

  # Argument 'cn':
  cn <- as.integer(cn);
  if (length(cn) == 1) {
    cn <- rep(cn, J);
  } else if (length(cn) != J) {
    stop("The length of argument 'cn' does not match 'y': ",
                                            length(cn), " != ", J);
  }
  uniqueCNs <- sort(unique(cn));
  unknown <- which(!is.element(uniqueCNs, c(0,1,2,NA)));
  if (length(unknown) > 0) {
    unknown <- paste(uniqueCNs[unknown], collapse=", ");
    stop("Argument 'cn' contains unknown CN levels: ", unknown);
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'adjust':
  adjust <- as.double(adjust);
  if (length(adjust) != 1) {
    stop("Argument 'adjust' must be single value: ", adjust);
  }
  if (adjust <= 0) {
    stop("Argument 'adjust' must be positive: ", adjust);
  }

  # Argument 'censorAt':
  censorAt <- Arguments$getDoubles(censorAt, length=c(2,2));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 
 

  verbose && enter(verbose, "Calling genotypes from allele B fractions (BAFs)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate result
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  naValue <- as.double(NA);
  mu <- rep(naValue, times=J);

  verbose && enter(verbose, "Censoring BAFs");
  verbose && cat(verbose, "Before:");
  verbose && summary(verbose, y);
  verbose && print(verbose, sum(is.finite(y)));
  # Censor values
  y[y < censorAt[1]] <- -Inf;
  y[y > censorAt[2]] <- +Inf;
  verbose && cat(verbose, "After:");
  verbose && summary(verbose, y);
  verbose && print(verbose, sum(is.finite(y)));
  verbose && exit(verbose);

  # To please R CMD check
  type <- NULL; rm(type);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call genotypes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq(along=uniqueCNs)) {
    cnKK <- uniqueCNs[kk];
    verbose && enter(verbose, sprintf("Copy number level #%d (C=%g) of %d", kk, cnKK, length(uniqueCNs)));

    keep <- which(cn == cnKK);
    yKK <- y[keep];

    # Exclude missing and non-finited values when fitting the density
    yT <- yKK[is.finite(yKK)];
    fit <- findPeaksAndValleys(yT, adjust=adjust, ...);
    verbose && cat(verbose, "Identified extreme points in density of BAF:");
    verbose && print(verbose, fit);

    fit <- subset(fit, type == "valley");
    nbrOfGenotypeGroups <- nrow(fit) + 1L;
    verbose && cat(verbose, "Local minimas (\"valleys\") in BAF:");
    verbose && print(verbose, fit);

    # Call genotypes
    muKK <- rep(naValue, length(yKK));
    if (cnKK == 0) {
      verbose && cat(verbose, "TCN=0 => BAF not defined. Skipping.");
    } else if (cnKK == 1) {
      verbose && cat(verbose, "TCN=1 => BAF in {0,1}.");
      # Sanity check
      stopifnot(nbrOfGenotypeGroups == 2);
      a <- fit$x[1];
      muKK[yKK <= a] <- 0;
      muKK[a < yKK] <- 1;
    } else if (cnKK == 2) {
      verbose && cat(verbose, "TCN=2 => BAF in {0,1/2,1}.");
      # Sanity check
      stopifnot(nbrOfGenotypeGroups == 3);
      a <- fit$x[1];
      b <- fit$x[2]; 
      muKK[yKK <= a] <- 0;
      muKK[a < yKK & yKK <= b] <- 1/2;
      muKK[b < yKK] <- 1;
    }
    mu[keep] <- muKK;

    verbose && exit(verbose);
  } # for (kk ...)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return genotype calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  stopifnot(length(mu) == J);

  verbose && exit(verbose);

  mu;
}) # callNaiveGenotypes()


###########################################################################
# HISTORY:
# 2009-11-03
# o Added an example() to the Rd help of callNaiveGenotypes().
# 2009-07-08
# o BUG FIX: Was never tested. Now tested via example(normalizeTumorBoost).
# 2009-07-06
# o Created from aroma.cn test script.
###########################################################################

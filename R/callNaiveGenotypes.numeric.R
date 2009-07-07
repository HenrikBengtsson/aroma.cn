###########################################################################/**
# @set "class=numeric"
# @RdocMethod callNaiveGenotypes
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
#  \item{...}{Not used.}
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
# @author
#*/########################################################################### 
setMethodS3("callNaiveGenotypes", "numeric", function(y, cn=rep(2L, length(y)), flavor=c("density"), adjust=1.5, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  J <- length(y);
  y <- as.double(y);

  # Argument 'cn':
  cn <- as.integer(cn);
  if (length(cn) != J) {
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate result
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  naValue <- as.double(NA);
  mu <- rep(naValue, times=J);


  # To please R CMD check
  type <- NULL; rm(type);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call genotypes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq(along=uniqueCNs)) {
    cnKK <- uniqueCNs[kk];
    keep <- which(cn == cnKK);
    yKK <- y[keep];
    muKK <- rep(naValue, length(yKK));
    fitKK <- findPeaksAndValleys(yKK, adjust=adjust);
    fit <- subset(fit, type == "valley");
    nbrOfGenotypeGroups <- nrow(fit)+1L;

    if (cnKK == 0) {
    } else if (cnKK == 1) {
      # Sanity check
      stopifnot(nbrOfGenotypeGroups == 2);
      a <- fit$x[1];
      muKK[yKK <= a] <- 0;
      muKK[a < yKK] <- 1;
    } else if (cnKK == 2) {
      # Sanity check
      stopifnot(nbrOfGenotypeGroups == 3);
      a <- fit$x[1];
      b <- fit$x[2]; 
      muKK[yKK <= a] <- 0;
      muKK[a < yKK & yKK <= b] <- 1/2;
      muKK[b < yKK] <- 1;
    }
    mu[keep] <- muKK;
  } # for (kk ...)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return genotype calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mu;
}) # callNaiveGenotypes()


###########################################################################
# HISTORY:
# 2009-07-06
# o Created from aroma.cn test script.
###########################################################################

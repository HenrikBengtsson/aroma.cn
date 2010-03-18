###########################################################################/**
# @set "class=numeric"
# @RdocMethod normalizeTumorBoost
# @alias normalizeTumorBoost
#
# @title "Normalizes allele B fractions for a tumor given a match normal"
#
# \description{
#  TumorBoost [1] is a normalization method that normalizes the allele B
#  fractions of a tumor sample given the allele B fractions and genotypes 
#  of a matched normal.
#  The method is a single-sample (single-pair) method. 
#  It does not require total copy-number estimates.
#  The normalization is done such that the total copy number is 
#  unchanged afterwards.
# }
# 
# @synopsis
#
# \arguments{
#  \item{betaT, betaN}{Two @numeric @vectors each of length J with 
#     tumor and normal allele B fractions, respectively.}
#  \item{muN}{An optional @vector of length J containing
#     normal genotypes calls in (0,1/2,1,@NA) for (AA,AB,BB).}
#  \item{flavor}{A @character string specifying the type of 
#     correction applied.}
#  \item{...}{Argument passed to @see "callNaiveGenotypes", if called.}
# }
#
# \value{
#   Returns a @numeric @vector of length J containing the normalized 
#   allele B fractions for the tumor.
#   Attribute \code{modelFit} is a @list containing model fit parameters.
# }
#
# \details{
#   Allele B fractions are defined as the ratio between the allele B signal
#   and the sum of both (all) allele signals at the same locus.
#   Allele B fractions are typically within [0,1], but may have a slightly 
#   wider support due to for instance negative noise.
#   This is typically also the case for the returned normalized 
#   allele B fractions.
# }
#
# \section{Flavors}{
#  This method provides a few different "flavors" for normalizing the
#  data.  The following values of argument \code{flavor} are accepted:
#  \itemize{
#   \item{v4: (default) The TumorBoost method, i.e. Eqns. (8)-(9) in [1].}
#   \item{v3: Eqn (9) in [1] is applied to both heterozygous and homozygous 
#             SNPs, which effectly is v4 where the normalized allele B
#             fractions for homozygous SNPs becomes 0 and 1.}
#   \item{v2: ...}
#   \item{v1: TumorBoost where correction factor is force to one, i.e.
#             \eqn{\eta_j=1}.  As explained in [1], this is a suboptimal
#             normalization method.  See also the discussion in the
#             paragraph following Eqn (12) in [1].}
#  }
# }
#
# @examples "../incl/normalizeTumorBoost.Rex"
#
# \author{Henrik Bengtsson and Pierre Neuvial}
#
# \references{
#  [1] H. Bengtsson, P. Neuvial & T.P. Speed, \emph{TumorBoost: Normalization
#      of allele-specific tumor copy numbers from a single pair of 
#      tumor-normal genotyping microarrays}, 2010 (submitted)\cr
# }
#*/########################################################################### 
setMethodS3("normalizeTumorBoost", "numeric", function(betaT, betaN, muN=callNaiveGenotypes(betaN), flavor=c("v4", "v3", "v2", "v1"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'betaT' & 'betaN':
  betaT <- as.numeric(betaT);
  betaN <- as.numeric(betaN);

  J <- length(betaT);
  if (length(betaN) != J) {
    stop("The length of arguments 'betaT' and 'betaN' differ: ", 
                                                   length(betaN), " != ", J);
  }
  
  # Argument: 'muN':
  if (length(muN) != J) {
    stop("Argument 'muN' does not match the number of loci: ", 
                                                   length(muN), " != ", J);
  }
  knownGenotypes <- c(0, 1/2, 1, NA);
  unknown <- which(!is.element(muN, knownGenotypes));
  n <- length(unknown);
  if (n > 0) {
    unknown <- unique(muN[unknown]);
    unknownStr <- paste(unknown, collapse=", ");
    stop("Argument 'muN' contains unknown values: ", unknownStr);
  }

  # Argument: 'flavor':
  flavor <- match.arg(flavor);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify set to be updated
  toUpdate <- which(is.finite(betaT) & is.finite(betaN) & is.finite(muN));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate delta
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  delta <- (betaN - muN);

  if (flavor == "v1") {
    b <- 1;
  } else if (flavor == "v2") {
    b <- rep(1, length(delta));
    isDown <- (betaT < betaN);
    isBetaNZero <- (betaN == 0);
    isBetaNOne <- (betaN == 1);
    idxs <- whichVector(isDown & !isBetaNZero);
    b[idxs] <- betaT[idxs]/betaN[idxs];
    idxs <- whichVector(!isDown & !isBetaNOne);
    b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
    rm(isDown, idxs);

    # Treat the case when the estimated SNP effect is zero
    # Then we want the normalized value to be exactly zero or one.
    idxs <- whichVector(delta == 0);
    
  } else if (flavor == "v3") {
    b <- rep(1, length(delta));
    isHomA <- (muN == 0);
    isHomB <- (muN == 1);
    isHet <- (!isHomA & !isHomB);
    isDown <- (betaT < betaN);
    isBetaNZero <- (betaN == 0);
    isBetaNOne <- (betaN == 1);
    idxs <- whichVector((isHet & isDown & !isBetaNZero) | (isHomA & !isBetaNZero));
    b[idxs] <- betaT[idxs]/betaN[idxs];
    idxs <- whichVector((isHet & !isDown & !isBetaNOne) | (isHomB & !isBetaNOne));
    b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
    rm(isDown, isHet, isHomA, isHomB, idxs);
  } else if (flavor == "v4") {
    # This is the published TumorBoost normalization method
    b <- rep(1, length(delta));
    isHet <- (muN != 0 & muN != 1);
    isDown <- (betaT < betaN);
    idxs <- whichVector(isHet & isDown);
    b[idxs] <- betaT[idxs]/betaN[idxs];
    idxs <- whichVector(isHet & !isDown);
    b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
    rm(isDown, isHet, idxs);
  }

  delta <- b * delta;

  # Sanity check
  stopifnot(length(delta) == J);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  betaTN <- betaT - delta;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return normalized data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  modelFit <- list(
    method = "normalizeTumorBoost",
    flavor = flavor,
    delta = delta
  );
  attr(betaTN, "modelFit") <- modelFit;

  # Sanity check
  stopifnot(length(betaTN) == J);

  betaTN;
}) # normalizeTumorBoost()


############################################################################
# HISTORY:
# 2010-03-18
# o BUG FIX: For flavors "v2" and "v3" NaN:s could be introduced if betaN
#   was exactly zero or exactly one.
# 2009-07-08
# o Now the arguments are 'betaT', 'betaN' and 'muN'.
# o Added an example() with real data.
# 2009-07-06
# o Created from process() of TumorBoostNormalization in aroma.cn.
# o Added model 'flavor' "v4" which corrects heterozygots according to "v2"
#   and homozygotes according to "v1".
# o Added model 'flavor' "v3".  Suggested by PN last night over a Guinness
#   at the pub after a long day of hard work.
# 2009-06-22
# o Added model 'flavor' "v2".
# 2009-06-08
# o The constructor of TumorBoostNormalization now only takes an
#   AromaUnitGenotypeCallSet for argument 'gcN'.  It no longer takes an
#   AromaUnitFracBCnBinarySet object.
# 2009-05-17
# o Now the constructor of TumorBoostNormalization asserts that there are
#   no stray arguments.
# 2009-04-29
# o Created.
############################################################################ 

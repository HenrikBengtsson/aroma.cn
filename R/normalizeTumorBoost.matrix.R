###########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizeTumorBoost
#
# @title "Normalizes allele B fractions for a tumor given a match normal"
#
# \description{
#  TumorBoost is normalization method that normalizes the allele B fractions
#  of a tumor sample given the allele B fractions and genotypes of a 
#  matched normal.
#  The method is a single-sample (single-pair) method. 
#  It does not require total copy-number estimates.
#  The normalization is done such that the total copy number is unchanged
#  afterwards.
# }
# 
# @synopsis
#
# \arguments{
#  \item{X}{A Jx2 @numeric @matrix of allele B fractions.  The first column
#     should contain values for the tumor and the second values for the
#     matched normal.}
#  \item{genotypes}{An optional @numeric @vector of length J containing
#     normal genotypes calls in (0,1/2,1) for (AA,AB,BB).}
#  \item{flavor}{A @character string specifying the type of 
#     correction applied.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric @vector of length J containing the normalized 
#   allele B fractions for the tumor.
#   Attribute \code{modelFit} is a @list containing model fit parameters.
# }
#
# \author{Henrik Bengtsson and Pierre Neuvial}
#*/########################################################################### 
setMethodS3("normalizeTumorBoost", "matrix", function(X, genotypes=NULL, flavor=c("v4", "v3", "v2", "v1"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'X':
  dim <- dim(X);
  if (dim[2] < 2) {
    stop("TumorBoost normalization requires two columns: ", dim[2]);
  }

  # Argument: 'genotypes':
  if (!is.null(genotypes)) {
    if (length(genotypes) != dim[1]) {
      stop("Argument 'genotypes' does not match the number of loci: ", 
                                          length(genotypes), " != ", dim[1]);
    }
    knownGenotypes <- c(0,1/2,1,NA);
    unknown <- which(!is.element(genotypes, knownGenotypes));
    n <- length(unknown);
    if (n > 0) {
      unknown <- unique(genotypes[unknown]);
      unknownStr <- paste(unknown, collapse=", ");
      stop("Argument 'genotypes' contains unknown values: ", unknownStr);
    }
  }

  # Argument: 'flavor':
  flavor <- match.arg(flavor);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  betaT <- X[,1];
  betaN <- X[,2];
  muN <- genotypes;
  
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
    idxs <- whichVector(isDown);
    b[idxs] <- betaT[idxs]/betaN[idxs];
    idxs <- whichVector(!isDown);
    b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
    rm(isDown,isHomA,isHomB,idxs);
  } else if (flavor == "v3") {
    b <- rep(1, length(delta));
    isHomA <- (muN == 0);
    isHomB <- (muN == 1);
    isHet <- !isHomA & !isHomB;
    isDown <- (betaT < betaN);
    idxs <- whichVector((isHet & isDown) | isHomA);
    b[idxs] <- betaT[idxs]/betaN[idxs];
    idxs <- whichVector((isHet & !isDown) | isHomB);
    b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
    rm(isDown,isHet,isHomA,isHomB,idxs);
  } else if (flavor == "v4") {
    b <- rep(1, length(delta));
    isHet <- (muN != 0 & muN != 1);
    isDown <- (betaT < betaN);
    idxs <- whichVector(isHet & isDown);
    b[idxs] <- betaT[idxs]/betaN[idxs];
    idxs <- whichVector(isHet & !isDown);
    b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
    rm(isDown,isHet,idxs);
  }
  delta <- b*delta;

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

  betaTN;
}) # normalizeTumorBoost()


############################################################################
# HISTORY:
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

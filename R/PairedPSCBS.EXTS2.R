###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod normalizeBAFsByRegions
#
# @title "Normalizes allele B fractions (BAFs) based on region-based PSCN estimates"
#
# \description{
#  @get "title" as given by the PSCBS segmentation method.
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A PairedPSCBS fit object as returned by 
#     @see "psCBS::segmentByPairedPSCBS".}
#   \item{by}{A @character string specifying if the normalization function
#     should be estimated based on TumorBoost normalized or non-normalized
#     tumor allele B fractions (BAFs).}
#   \item{...}{Additional arguments passed
#     @see "aroma.cn::normalizeMirroredBAFsByRegions".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a PairedPSCBS fit object where the region-level 
#   decrease-in-heterozygosity (DH) means have been normalized,
#   as well as the locus-specific tumor allele B fractions.
# }
#
# \details{
#   Note that his normalization method depends on the segmentation
#   results. Hence, it recommended \emph{not} to resegment the
#   normalized signals returned by this, because such a segmentation
#   will be highly dependent on the initial segmentation round.
# }
#
# @examples "../incl/normalizeBAFsByRegions.PairedPSCBS.Rex"
#
# @author
#
# \seealso{
#   Internally @see "aroma.cn::normalizeMirroredBAFsByRegions" is used.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("normalizeBAFsByRegions", "PairedPSCBS", function(fit, by=c("betaTN", "betaT"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'by':
  by <- match.arg(by);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Normalizes region-level mirrored allele B fractions (mBAFs)");

  data <- fit$data;
  stopifnot(!is.null(data));

  segs <- fit$output;
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  x <- data$x;
  betaT <- data$betaT;
  betaTN <- data$betaTN;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate region-level mBAFs for homozygous SNPs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calculating region-level mBAFs for homozygous SNPs");

  # Calculate mBAFs for all loci
  beta <- data[[by]];
  rho <- 2*abs(beta - 1/2);
  rm(beta);

  # Identify homozygous SNPs
  muN <- data$muN;
  isHom <- (muN == 0 | muN == 1);
  rm(muN);

  # Allocate
  naValue <- as.double(NA);
  X <- matrix(naValue, nrow=nbrOfSegments, ncol=3);
  for (kk in seq(length=nbrOfSegments)) {
    xRange <- as.numeric(segs[kk,c("dh.loc.start", "dh.loc.end")]);
    tcn <- segs[kk,"tcn.mean"];
    dh <- segs[kk,"dh.mean"];
    # Identify all homozygous SNPs in the region
    idxs <- whichVector(xRange[1] <= x & x <= xRange[2] & isHom);
    mBAFhom <- mean(rho[idxs], na.rm=TRUE);
    X[kk,] <- c(dh, mBAFhom, tcn);
  } # for (kk ...)

  # Not needed anymore
  rm(rho, isHom);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Normalize region-level DHs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Normalizing region-level mBAFs");
  XN <- normalizeMirroredBAFsByRegions(data=X, ..., verbose=verbose);

  # Update DH segmentation means
  rhoN <- XN[,1,drop=TRUE];
  segs[,"dh.mean"] <- rhoN;
  rm(rhoN);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Normalize locus-level data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Normalizing locus-level data accordingly");
  modelFit <- attr(XN, "modelFit");
  scale <- modelFit$scale;
  rm(modelFit, XN);

  # Expand region-level scale factors to locus-level scale factors
  naValue <- as.double(NA);
  scales <- rep(naValue, times=length(data$betaT));
  for (kk in seq(length=nbrOfSegments)) {
    xRange <- as.numeric(segs[kk,c("dh.loc.start", "dh.loc.end")]);
    # Identify all SNPs in the region
    idxs <- whichVector(xRange[1] <= x & x <= xRange[2]);
    scales[idxs] <- scale[kk];
  } # for (kk ...)

  # Update tumor allele B fractions
  for (ff in c("betaT", "betaTN")) {
    beta <- data[[ff]];
    beta <- beta - 1/2;
    beta <- scales * beta;
    beta <- beta + 1/2;
    data[[ff]] <- beta;
  }
  rm(beta, scales);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fitN <- fit;
  fitN$data <- data;
  fitN$output <- segs;

  verbose && exit(verbose);

  fitN;
}) # normalizeBAFsByRegions()


###########################################################################/**
# @RdocMethod callAllelicBalanceByBAFs
#
# @title "Calls regions that are in allelic balance"
#
# \description{
#  @get "title" from the allele B fractions (BAF).
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A PairedPSCBS fit object as returned by 
#     @see "psCBS::segmentByPairedPSCBS".}
#   \item{maxScore}{A positive @double threshold.}
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, an already called object is skipped, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a PairedPSCBS fit object where columns for
#   allelic imbalance scores and p-values as well as allelic
#   balance calls are added.
# }
#
# @examples "../incl/callAllelicBalanceByBAFs.PairedPSCBS.Rex"
#
# @author
#
# @keyword internal
#*/########################################################################### 
setMethodS3("callAllelicBalanceByBAFs", "PairedPSCBS", function(fit, maxScore=4, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'maxScore':
  maxScore <- Arguments$getDouble(maxScore, range=c(0,Inf));

  # Extract segments
  segs <- as.data.frame(fit);


  # Nothing to do?
  if (!force && !is.null(segs$ab.call)) {
    # Allelic balance segments are already called
    return(fit);
  }


  # Extract data
  betaTN <- fit$data$betaTN;
  muN <- fit$data$muN;

  nbrOfSegments <- nrow(segs);
  naValue <- as.double(NA);
  df <- NULL;
  for (kk in seq(length=nbrOfSegments)) {
    fitS <- subsetBySegments(fit, idxs=kk);
    dataS <- fitS$data;
    betaTN <- dataS$betaTN;
    muN <- dataS$muN;
    fitKK <- testAllelicBalanceByBAFs(betaTN, muN=muN);

    dfKK <- data.frame(
      statistic=fitKK$statistic,
      p.value=fitKK$p.value
    );

    df <- rbind(df, dfKK);
  } # for (kk ...)
  rownames(df) <- NULL;

  colnames(df) <- c("ai", "ai.p.value");
  df$ab.call <- (df$ai <= maxScore);

  segs <- cbind(segs, df);

  fitC <- fit;
  fitC$output <- segs;

  fitC;
})



###########################################################################/**
# @RdocMethod callCopyNeutralRegions
#
# @title "Calls regions that are copy neutral"
#
# \description{
#  @get "title" from the allele B fractions (BAF).
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A PairedPSCBS fit object as returned by 
#     @see "psCBS::segmentByPairedPSCBS".}
#   \item{...}{Additional arguments passed to 
#     @see "aroma.cn::findNeutralCopyNumberState".}
#   \item{force}{If @TRUE, an already called object is skipped, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a PairedPSCBS fit object where a column with the copy-neutral call.
# }
#
# @examples "../incl/callCopyNeutralRegions.PairedPSCBS.Rex"
#
# @author
#
# @keyword internal
#*/########################################################################### 
setMethodS3("callCopyNeutralRegions", "PairedPSCBS", function(fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  # Call allelic balance or not, unless already done
  fit <- callAllelicBalanceByBAFs(fit, ..., verbose=verbose);

  segs <- fit$output;

  # Nothing to do?
  if (!force && !is.null(segs$neutral.call)) {
    # Copy neutral segments are already called
    return(fit);
  }

  C <- segs[,"tcn.mean", drop=TRUE];
  isAB <- segs[,"ab.call", drop=TRUE];
  n <- segs[,"tcn.num.snps", drop=TRUE]; # "tcn.num.mark"? /HB 2010-09-09

  # Give more weight to longer regions
  weights <- n;

  isNeutral <- findNeutralCopyNumberState(C=C, isAI=!isAB, weights=weights,
                                                       ..., verbose=verbose);

  segs$neutral.call <- isNeutral;

  fitC <- fit;
  fitC$output <- segs;
  
  fitC;
})




##############################################################################
# HISTORY
# 2010-09-15 [HB]
# o Added Rdocs for callCopyNeutralRegions().
# 2010-09-09 [HB]
# o Added callCopyNeutralRegions() for PairedPSCBS.
# 2010-09-08 [HB]
# o Added subsetBySegments() for PairedPSCBS.
# o Added Rdocs with an example.
# o Added normalizeBAFsByRegions() for PairedPCSBS.
# o Created.
##############################################################################

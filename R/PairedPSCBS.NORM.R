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
# @set "class=PairedPSCBS"
# @RdocMethod orthogonalizeC1C2
#
# @title "Orthogonalizes (C1,C2) based on region-based PSCN estimates"
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
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a PairedPSCBS fit object.
# }
#
# @examples "../incl/orthogonalizeC1C2.PairedPSCBS.Rex"
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("orthogonalizeC1C2", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Orthogonalizes region-level minor and major CNs (C1,C2)");

  data <- fit$data;
  stopifnot(!is.null(data));

  segs <- as.data.frame(fit);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (C1,C2,...)
  X <- extractMinorMajorCNs(fit);

  # Number of TCN and DH data points
  counts <- X[,3:4, drop=FALSE];

  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Orthogonalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (dC1,dC2)
  dC1C2 <- colDiffs(C1C2);

  # Slopes
  beta <- atan(dC1C2[,1]/dC1C2[,2]);

  # Find mode in {beta} closest to zero
  fp <- findPeaksAndValleys(beta);
  type <- NULL; rm(type); # To please R CMD check
  pfp <- subset(fp, type == "peak");
  ww <- which.min(abs(pfp[,"x"])); # closest to 0
  bb <- pfp[ww,"x"];

  # Verticalization
  C1 <- C1C2[,1];
  C2 <- C1C2[,2];
  C1 <- C1-sin(bb)*C2;  ## projection...
  C1C2o <- C1C2;
  C1C2o[,1] <- C1;
  C1C2o[,2] <- C2;

  # (C1,C2) -> (TCN,DH)
  gamma <- C1 + C2;
  dh <- C2 / gamma;

  # Update segmentation means
  segs[,"tcn.mean"] <- gamma;
  segs[,"dh.mean"] <- dh;

  # Update data [TO DO]


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fitO <- fit;
  fitO$output <- segs;

  verbose && exit(verbose);

  fitO;
}) # orthogonalizeC1C2()



##############################################################################
# HISTORY
# 2010-09-19 [HB+PN]
# o Added orthogonalizeC1C2() for PairedPSCBS.
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

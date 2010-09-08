###########################################################################/**
# @RdocDefault testAllelicImbalanceByAlleleBFractions
#
# @title "Tests for allelic imbalance in a genomic region"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{betaT}{A @numeric @vector of J tumor allele B fractions.}
#   \item{muN}{A @numeric @vector of J normal (germline) genotypes.}
#   \item{flavor}{A @character specifying the type of test to be performed.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   A @list of class "htest".
# }
#
# @author
#*/########################################################################### 
setMethodS3("testAllelicImbalanceByAlleleBFractions", "default", function(betaT, muN, flavor=c("var", "mean"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'betaT':
  betaT <- Arguments$getDoubles(betaT);
  nbrOfLoci <- length(betaT);
  length2 <- rep(nbrOfLoci, times=2);

  # Argument 'muN':
  muN <- Arguments$getDoubles(muN, length=length2);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 
 

  verbose && enter(verbose, "Testing for allelic imbalance by allele B fractions");

  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);
  verbose && cat(verbose, "Requested flavor: ", flavor);
  verbose && cat(verbose, "Tumor allele B fractions:");
  verbose && str(verbose, betaT);
  verbose && cat(verbose, "Match-normal genotypes:");
  verbose && str(verbose, muN);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop loci with missing signals or missing genotypes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Dropping loci with missing signals or missing genotypes");
  keep <- (is.finite(betaT) & is.finite(muN));
  keep <- whichVector(keep);

  if (length(keep) < nbrOfLoci) {
    betaT <- betaT[keep];
    muN <- muN[keep];
    nbrOfLoci <- length(betaT);
    verbose && cat(verbose, "Number of SNPs with finite signals: ", nbrOfLoci);
    verbose && cat(verbose, "Tumor allele B fractions:");
    verbose && str(verbose, betaT);
    verbose && cat(verbose, "Match-normal genotypes:");
    verbose && str(verbose, muN);
  }

  # Not needed anymore
  rm(keep);
  verbose && exit(verbose);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get homozygous and heterozygous centered BAF signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Getting homozygous and heterozygous centered BAF signals");
  # AA SNPs
  isAA <- (muN == 0);
  idxs <- whichVector(isAA);
  yAA <- betaT[idxs];
  yAA <- yAA - mean(yAA);
  
  # BB SNPs
  isBB <- (muN == 1);
  idxs <- whichVector(isBB);
  yBB <- betaT[idxs];
  yBB <- yBB - mean(yBB);

  # AB SNPs
  isAB <- (muN == 1/2);
  idxs <- whichVector(isAB);
  yAB <- betaT[idxs];
  yAB <- yAB - mean(yAB);

  # Homozygous SNPs
  homs <- c(yAA, yBB);
  nHoms <- length(homs);

  # Heterozygous SNPs
  hets <- yAB;
  nHets <- length(hets);

  # All other SNPs (for curiousity only)
  ## isOther <- (!isAA & !isAB & !isBB);
  ## idxs <- whichVector(isOther);
  ## nOthers <- length(idxs);

  # Not needed anymore
  rm(isAA, isBB, isAB, idxs, yAA, yBB, yAB);

  verbose && cat(verbose, "Number of homozygous SNPs: ", nHoms);
  verbose && cat(verbose, "Number of heterozygous SNPs: ", nHets);

  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nMin <- min(nHoms, nHets);

  # Sanity check
  if (nMin <= 1) {
    throw(sprintf("Cannot infer allelic imbalance. Too few (homozygous, heterozygous) data points: (%d,%d)", nHoms, nHets));
  }

  # Need to fall back on another test than the requested one?
  if (flavor == "var") {
    if (nMin <= 2) {
      warning(sprintf("Cannot infer allelic imbalance using the heterozygous-homozygous variance test. Will use the mean test instead. Too few (homozygous, heterozygous) data points: (%d,%d)", nHoms, nHets));
      flavor <- "mean";
     verbose && cat(verbose, "Adjusted flavor: ", flavor);
    }
  } else if (flavor == "mean") {
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Testing for allelic imbalance
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Default statistics
  res <- list(statistic=NA, p.value=NA);

  if (flavor == "var") {
    res <- var.test(hets, homs, alternative="greater");
  } else if (flavor == "mean") {
    sigmaHat <- sd(homs, na.rm=TRUE)
    stat <- mean(abs(hets), na.rm=TRUE) / sigmaHat;
    p <- 1 - pnorm(stat);
    res <- list(statistic=stat, p.value=p);
    class(res) <- c("htest", class(res));
  }

  verbose && str(verbose, res);

  verbose && exit(verbose);
  
  res;
}, protected=TRUE) # testAllelicImbalanceByAlleleBFractions()


#############################################################################
# HISTORY
# 2010-09-08 [PN,HB]
# o Harmonized code. Added more verbose output. Added more code comments.
# o Renamed to testAllelicImbalanceByAlleleBFractions() from AI.test().
# o Added an Rdoc skeleton.
# o Renamed arguments. 
# 2010-08-25 [PN]
# o Created.
#############################################################################

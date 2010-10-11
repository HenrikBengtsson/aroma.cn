###########################################################################/**
# @set "class=PairedPSCBS"
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
setMethodS3("callAllelicBalanceByBAFs", "PairedPSCBS", function(fit, maxScore=4, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'maxScore':
  maxScore <- Arguments$getDouble(maxScore, range=c(0,Inf));

  # Extract segments
  segs <- as.data.frame(fit);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'cache':
  cache <- Arguments$getLogical(cache);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Nothing to do?
  if (!force && !is.null(segs$ab.call)) {
    # Allelic balance segments are already called
    return(fit);
  }

  verbose && enter(verbose, "Calling allelic balance by BAFs");

  # Extract data
  betaTN <- fit$data$betaTN;
  muN <- fit$data$muN;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="callAllelicBalanceByBAFs", class=class(fit)[1], 
    data=list(betaTN=betaTN, muN=muN, segments=as.data.frame(fit)),
    maxScore = maxScore,
    version="2010-10-10"
  );
  dirs <- c("aroma.cn", "ortho");
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(res);
    }
  }



  nbrOfSegments <- nrow(segs);
  naValue <- as.double(NA);
  df <- NULL;
  for (kk in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment #%d of %d", kk, nbrOfSegments));

    fitS <- subsetBySegments(fit, idxs=kk);
    verbose && print(verbose, fitS);

    dataS <- fitS$data;
    betaTN <- dataS$betaTN;
    muN <- dataS$muN;

    verbose && summary(verbose, betaTN);
    verbose && summary(verbose, muN);

    # AD HOC: For some unknown reason does resample() introduce NAs.
    # /HB 2010-09-15
    keep <- is.finite(betaTN) & is.finite(muN);
    keep <- whichVector(keep);
    betaTN <- betaTN[keep];
    muN <- muN[keep];

    fitKK <- testAllelicBalanceByBAFs(betaTN, muN=muN);

    dfKK <- data.frame(
      statistic=fitKK$statistic,
      p.value=fitKK$p.value
    );

    df <- rbind(df, dfKK);

    verbose && exit(verbose);
  } # for (kk ...)
  rownames(df) <- NULL;

  colnames(df) <- c("ai", "ai.p.value");
  df$ab.call <- (df$ai <= maxScore);

  segs <- cbind(segs, df);

  fitC <- fit;
  fitC$output <- segs;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    saveCache(key=key, dirs=dirs, fitC);
  }

  verbose && exit(verbose);

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


setMethodS3("plotC1C2Grid", "PairedPSCBS", function(fit, ..., Clim=c(0,4), main=NULL) {
  plotC1C2(fit, Clim=Clim);
  title(main=main);
  data <- extractC1C2(fit);
  n <- data[,4, drop=TRUE];
  n <- sqrt(n);
  w <- n/sum(n, na.rm=TRUE);
  adjust <- 0.2;
  for (cc in 1:2) {
    y <- data[,cc];
    ok <- is.finite(y) & is.finite(w);
    y <- y[ok];
    wt <- w[ok]/sum(w[ok]);
    d <- density(y, weights=wt, adjust=adjust);
    draw(d, side=cc, height=0.3, col="gray", lwd=2, xpd=FALSE);
    if (cc == 2) {
      draw(d, side=1, height=0.3, col="lightblue", lwd=2, xpd=FALSE);
    }
    p <- findPeaksAndValleys(d, tol=0.05);
    type <- NULL; rm(type); # To please R CMD check
    p <- subset(p, type == "peak");
    p <- p[order(p$density, decreasing=TRUE),,drop=FALSE];
    p <- head(p, n=8);
    if (cc == 1) {
      abline(v=p$x, lty=3, col="gray");
    } else {
      abline(h=p$x, lty=3, col="gray");
    }
  }
  box();
})

##############################################################################
# HISTORY
# 2010-10-10 [HB]
# o Added memoization to callAllelicBalanceByBAFs().
# 2010-10-08 [HB]
# o Added plotC1C2Grid() for PairedPSCBS.
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

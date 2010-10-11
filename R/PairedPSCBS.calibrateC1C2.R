###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod calibrateC1C2
#
# @title "Calibrates ASCN signals in (C1,C2) space based on region-based PSCN estimates"
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
#   Returns a calibrated PairedPSCBS fit object.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("calibrateC1C2", "PairedPSCBS", function(fit, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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


  verbose && enter(verbose, "Calibrating ASCNs in (C1,C2) space");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update TCN segment means
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  key <- list(method="postsegmentTCN", class=class(fit)[1], 
    data=as.data.frame(fit),
    version="2010-10-10"
  );
  dirs <- c("aroma.cn", "ortho");
  if (!force) {
    fit2 <- loadCache(key=key, dirs=dirs);
    if (!is.null(fit2)) {
      verbose && cat(verbose, "Cached results found.");
    }
  }

  if (is.null(fit2)) {
    fit2 <- postsegmentTCN(fit, verbose=verbose);
    if (cache) {
      saveCache(key=key, dirs=dirs, fit2);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove effects in BAF that are dependent on TCN.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit3 <- normalizeBAFsByRegions(fit2, force=force, cache=cache, verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call regions in allelic balance and shift them to (C1 = C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit4 <- callAllelicBalanceByBAFs(fit3, verbose=verbose);
  ww <- which(fit4$output$ab.call);
  fit4$output[ww, "dh.mean"] <- 0;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Deshear by (C1,C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit5 <- deShearC1C2(fit4, dirs="|-", verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Deshear by (C1,C2) - diagonals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit6 <- deShearC1C2(fit5, dirs="X", verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Deshear by (C1,C2) - vertical, horizontal, diagonals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit7 <- deShearC1C2(fit6, dirs="|,-,X");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove offset in C2
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dC2 <- estimateC2Bias(fit7, verbose=verbose);
  # Sanity check
  dC2 <- Arguments$getDouble(dC2, range=c(-3,3));
  fit8 <- translateC1C2(fit7, dC2=-dC2, verbose=verbose);


  verbose && exit(verbose);

 
  fit8;
})



##############################################################################
# HISTORY
# 2010-10-10 [HB]
# o Added calibrateC1C2() for PairedPSCBS.
##############################################################################

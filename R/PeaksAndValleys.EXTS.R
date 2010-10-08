setMethodS3("callPeaks", "data.frame", function(fit, ...) {
  # Argument 'fit';
  stopifnot(all(is.element(c("type", "x", "density"), colnames(fit))));
  class(fit) <- c("PeaksAndValleys", class(fit));
  callPeaks(fit, ...);
}) # callPeaks()


setMethodS3("callPeaks", "PeaksAndValleys", function(fit, expected=c(-1/2,-1/4,0,+1/4,+1/2)*pi, ..., flavor=c("decreasing", "all"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fit';
  stopifnot(all(is.element(c("type", "x", "density"), colnames(fit))));

  # Argument 'expected':
  expected <- Arguments$getNumerics(expected);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling peaks");
  verbose && cat(verbose, "Flavor: ", flavor);

  verbose && cat(verbose, "All expected peaks:");
  verbose && print(verbose, expected);

  verbose && enter(verbose, "Extracing peaks");
  subset <- which(fit$type == "peak");
  fitP <- fit[subset,,drop=FALSE];
  verbose && print(verbose, fitP);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calling peaks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  xd <- fitP[,c("x", "density"),drop=FALSE];

  if (flavor == "all") {
    calls <- sapply(xd$x, FUN=function(x) {
      dist <- abs(x - expected);
      which.min(dist);
    });
  } else if (flavor == "decreasing") {
    # It is probably better to call the strongest peaks first for which
    # we have more confidence, and then call the other relative to those.
    # /HB 2010-09-19

    # Order peaks by density
    o <- order(xd[,"density"], decreasing=TRUE);
    verbose && cat(verbose, "Reordering:");
    verbose && print(verbose, o);
    # The ranks (for later)
    r <- seq(along=o); r[o] <- r;
    verbose && cat(verbose, "Rank:");
    verbose && print(verbose, r);
    xd <- xd[o,,drop=FALSE];
    verbose && print(verbose, xd);

    # Call the strongest peak first, then the 2nd strongest and so on...
    naValue <- as.integer(NA);
    calls <- rep(naValue, times=nrow(xd));
    expectedLeft <- expected;
    for (kk in seq(length=nrow(xd))) {
      # All expected modes called?
      if (!any(is.finite(expectedLeft))) {
        break;
      }
      # Mode #kk
      x <- xd[kk,"x"];
      dx <- abs(x - expectedLeft);
      call <- which.min(dx);
      expectedLeft[call] <- NA;
      calls[kk] <- call;
    } # for (kk ...)
  } # if (flavor ...)

  verbose && cat(verbose, "Calls:");
  verbose && print(verbose, calls);
  verbose && cat(verbose, "Expected values:");
  verbose && print(verbose, expected[calls]);

  fitC <- cbind(fit, callId=as.integer(NA), call=as.double(NA));
  fitC[subset,"callId"] <- calls[r];
  fitC[subset,"call"] <- expected[calls[r]];
  attr(fitC, "expected") <- expected;

  verbose && print(verbose, fitC);

  verbose && exit(verbose);

  fitC;
}, protected=TRUE) # callPeaks()


##############################################################################
# HISTORY
# 2010-10-08 [HB]
# o Added callPeaks().
# o Created.
##############################################################################

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
# \examples{\dontrun{
#  @include "../incl/calibrateC1C2.PairedPSCBS.Rex"
# }}
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Align C2 peaks to C1 peaks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  peakFit <- fitC1C2Peaks(fit8);
  a <- peakFit$params$a;
  b <- peakFit$params$b;
  scale <- 1/b;
  shift <- -a/b;

  # Sanity check
  shift <- Arguments$getDouble(shift, range=c(-3,3));
  scale <- Arguments$getDouble(scale, range=c(0.1,3));
  fit9 <- translateC1C2(fit8, sC2=scale, dC2=shift, verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove offset in (C1,C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Removing offset (C1,C2) space");
  ff <- fitC1C2Densities(fit9);
  offset <- ff$pList$C1$x[1];
  verbose && cat(verbose, "Offset: ", offset);

  shift <- -offset;
  # Sanity check
  shift <- Arguments$getDouble(shift, range=c(-3,3));

  fit10 <- translateC1C2(fit9, dC1=shift, dC2=shift, verbose=verbose);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rescale (1,1) in (C1,C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Removing offset (C1,C2) space");
  ff <- fitC1C2Densities(fit10);
  scale <- ff$pList$C1$x[2];
  verbose && cat(verbose, "Scale: ", scale);

  scale <- 1/scale;
  # Sanity check
  scale <- Arguments$getDouble(scale, range=c(0.1,10));

  fit11 <- translateC1C2(fit10, sC1=scale, sC2=scale, verbose=verbose);
  verbose && exit(verbose);


  verbose && exit(verbose);

 
  fit11;
})



setMethodS3("fitC1C2Peaks", "PairedPSCBS", function(fit, ..., tol=0.05) {
  d1d2 <- fitC1C2Densities(fit);
  pList <- d1d2$pList;
  D <- outer(pList$C1$x, pList$C2$x, FUN="-");
  idxs <- which(abs(D) < tol, arr.ind=TRUE);

  if (length(idxs) < 2) {
    throw("Cannot fit relationship. Too few common peaks: ", length(idxs));
  }

  c1 <- pList$C1$x[idxs[,1]];
  c2 <- pList$C2$x[idxs[,2]];
  dd <- cbind(pList$C1$density[idxs[,1]], pList$C2$density[idxs[,2]]);
  dd <- rowSums(dd^2);
  w <- dd / sum(dd);
  f <- lm(c2 ~ 1 + c1, weights=w);
  print(f);
  a <- coef(f)[1];
  b <- coef(f)[2];
  params <- list(a=a, b=b);
  list(fit=f, params=params);
}) # fitC1C2Peaks()


setMethodS3("fitC1C2Densities", "PairedPSCBS", function(fit, adjust=0.2, tol=0.05, ...) {
  data <- extractC1C2(fit);
  n <- data[,4, drop=TRUE];
  n <- sqrt(n);
  w <- n/sum(n, na.rm=TRUE);
  adjust <- 0.2;

  dList <- list();
  for (cc in 1:2) {
    y <- data[,cc];
    ok <- is.finite(y) & is.finite(w);
    y <- y[ok];
    wt <- w[ok]/sum(w[ok]);
    d <- density(y, weights=wt, adjust=adjust);
    dList[[cc]] <- d;
  }
  names(dList) <- colnames(data)[1:2];

  type <- NULL; rm(type);  # To please R CMD check
  pList <- lapply(dList, FUN=function(d) {
    p <- findPeaksAndValleys(d, tol=tol);
    p <- subset(p, type == "peak");
    p <- p[order(p$density, decreasing=TRUE),,drop=FALSE];
  });
  names(pList) <- names(dList);

  return(list(dList=dList, pList=pList));
}) # fitC1C2Densities() 


##############################################################################
# HISTORY
# 2010-10-10 [HB]
# o Added fitC1C2Peaks().
# o Added calibrateC1C2() for PairedPSCBS.
##############################################################################

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
#     @see "PSCBS::segmentByPairedPSCBS".}
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
setMethodS3("calibrateC1C2", "PairedPSCBS", function(fit, ..., force=FALSE, cache=TRUE, debug=FALSE, verbose=FALSE) {
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
  } else {
    fit2 <- NULL;
  }

  if (is.null(fit2)) {
    verbose && enter(verbose, "Updating TCN statistics per DH segment");
    fit2 <- postsegmentTCN(fit, verbose=verbose);
    if (cache) {
      saveCache(key=key, dirs=dirs, fit2);
    }
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove effects in BAF that are dependent on TCN.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing BAFs by segment");
  fit3 <- normalizeBAFsByRegions(fit2, force=force, cache=cache, verbose=verbose);
  verbose && exit(verbose);

  if (debug) {
    ff <- fit3;
    figName <- "debug,fit3";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call regions in allelic balance and shift them to (C1 = C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Adjusting for biases in allelic-balance segments");
  fit4 <- callAllelicBalanceByBAFs(fit3, verbose=verbose);
  segs4 <- getSegments(fit4, splitters=TRUE);
  ww <- which(segs4$abCall);
  fit4$output[ww, "dhMean"] <- 0;
  verbose && exit(verbose);

  if (debug) {
    ff <- fit4;
    figName <- "debug,fit4";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Deshear by (C1,C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit5 <- deShearC1C2(fit4, dirs="|-", verbose=verbose);
  if (debug) {
    ff <- fit5;
    figName <- "debug,fit5";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Deshear by (C1,C2) - diagonals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit6 <- deShearC1C2(fit5, dirs="X", verbose=verbose);
  if (debug) {
    ff <- fit6;
    figName <- "debug,fit6";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Deshear by (C1,C2) - vertical, horizontal, diagonals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit7 <- deShearC1C2(fit6, dirs="|,-,X", verbose=verbose);
  if (debug) {
    ff <- fit7;
    figName <- "debug,fit7";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove offset in C2
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dC2 <- estimateC2Bias(fit7, verbose=verbose);
  # Sanity check
  dC2 <- Arguments$getDouble(dC2, range=c(-3,3));
  fit8 <- translateC1C2(fit7, dC2=-dC2, verbose=verbose);
  if (debug) {
    ff <- fit8;
    figName <- "debug,fit8";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Align C2 peaks to C1 peaks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  peakFit <- fitC1C2Peaks(fit8, onError="skip", verbose=verbose);
  if (!is.null(peakFit)) {
    a <- peakFit$params$a;
    b <- peakFit$params$b;
    scale <- 1/b;
    shift <- -a/b;
  
    # Sanity check
    shift <- Arguments$getDouble(shift, range=c(-3,3));
    scale <- Arguments$getDouble(scale, range=c(0.1,3));
    fit9 <- translateC1C2(fit8, sC2=scale, dC2=shift, verbose=verbose);
  } else {
    fit9 <- fit8;
  }

  if (debug) {
    ff <- fit9;
    figName <- "debug,fit9";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call regions in allelic balance and shift them to (C1 = C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Adjusting for biases in allelic-balance segments");
  fit10 <- callAllelicBalanceByBAFs(fit9, force=TRUE, verbose=verbose);
  segs10 <- getSegments(fit10, splitters=TRUE);
  ww <- which(segs10$abCall);
  fit10$output[ww, "dhMean"] <- 0;
  verbose && exit(verbose);

  if (debug) {
    ff <- fit10;
    figName <- "debug,fit9b";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }

  fit9 <- fit10;



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove offset in (C1,C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Removing offset (C1,C2) space");
  ff <- fitC1C2Densities(fit9, orderBy="x");
  pp <- ff$pList$C1;
  pp <- subset(pp, density > 0.5);
  offset <- sort(pp$x)[1];
  verbose && cat(verbose, "Offset: ", offset);

  shift <- -offset;
  # Sanity check
  shift <- Arguments$getDouble(shift, range=c(-3,3));

  fit10 <- translateC1C2(fit9, dC1=shift, dC2=shift, verbose=verbose);

  if (debug) {
    ff <- fit10;
    figName <- "debug,fit10";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rescale (1,1) in (C1,C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Rescaling (1,1) cluster in (C1,C2) space");
  ff <- fitC1C2Densities(fit10, orderBy="x");
  pp <- ff$pList$C1;
  pp <- subset(pp, density > 0.5);
  scale <- sort(pp$x)[2];
  verbose && cat(verbose, "Scale: ", scale);

  scale <- 1/scale;
  # Sanity check
  scale <- Arguments$getDouble(scale, range=c(0.1,10));

  fit11 <- translateC1C2(fit10, sC1=scale, sC2=scale, verbose=verbose);

  if (debug) {
    ff <- fit11;
    figName <- "debug,fit11";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Align C2 peaks to C1 peaks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  peakFit <- fitC1C2Peaks(fit11, tol=0.1, onError="skip", verbose=verbose);
  if (!is.null(peakFit)) {
    a <- peakFit$params$a;
    b <- peakFit$params$b;
    scale <- 1/b;
    shift <- -a/b;
  
    # Sanity check
    shift <- Arguments$getDouble(shift, range=c(-3,3));
    scale <- Arguments$getDouble(scale, range=c(0.1,3));
    fit12 <- translateC1C2(fit11, sC2=scale, dC2=shift, verbose=verbose);
  } else {
    fit12 <- fit11;
  }

  if (debug) {
    ff <- fit12;
    figName <- "debug,fit12";
    devSet(figName); devSet(figName);
    plotC1C2Grid(ff); linesC1C2(ff); stext(side=3,pos=1,figName);
  }


  verbose && exit(verbose);

  fit12;
})



setMethodS3("fitC1C2Peaks", "PairedPSCBS", function(fit, ..., tol=0.05, onError=c("error", "warning", "skip"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'onError':
  onError <- match.arg(onError);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting the relationship between peaks in C1 and C2");

  d1d2 <- fitC1C2Densities(fit, orderBy="x");
  pList <- d1d2$pList;
  verbose && print(verbose, pList);

  D <- outer(pList$C1$x, pList$C2$x, FUN="-");
  idxs <- which(abs(D) < tol, arr.ind=TRUE);
  verbose && print(verbose, idxs);

  if (nrow(idxs) < 2) {
    msg <- sprintf("Cannot fit relationship between C1 and C2. Too few common peaks: %d", nrow(idxs));
    if (onError == "error") {
      devSet("Exception"); devSet("Exception"); plotC1C2Grid(fit);
      throw(msg);
    }
    if (onError == "warning") {
      warning(msg);
    }

    msg <- sprintf("Skipping fitting the relationship between peaks in C1 and C2. %s", msg);
    verbose && cat(verbose, msg);
    warning(msg);
    verbose && exit(verbose);
    return(NULL);
  }

  c1 <- pList$C1$x[idxs[,1]];
  c2 <- pList$C2$x[idxs[,2]];
  dd <- cbind(pList$C1$density[idxs[,1]], pList$C2$density[idxs[,2]]);
  dd <- rowSums(dd^2);
  w <- dd / sum(dd);
  verbose && print(verbose, cbind(c1=c1, c2=c2, weights=w));

  f <- lm(c2 ~ 1 + c1, weights=w);
  verbose && print(verbose, f);

  a <- coef(f)[1];
  b <- coef(f)[2];
  params <- list(a=a, b=b);
  verbose && print(verbose, params);
  res <- list(fit=f, params=params);

  verbose && str(verbose, res);

  verbose && exit(verbose);

  res;
}) # fitC1C2Peaks()


setMethodS3("fitC1C2Densities", "PairedPSCBS", function(fit, adjust=0.2, tol=0.05, orderBy=c("density", "x"), ...) {
  orderBy <- match.arg(orderBy);

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
    p <- p[order(p[[orderBy]], decreasing=c("x"=FALSE, "density"=TRUE)[orderBy]),,drop=FALSE];
  });
  names(pList) <- names(dList);

  return(list(dList=dList, pList=pList));
}) # fitC1C2Densities() 


##############################################################################
# HISTORY
# 2011-10-16 [HB]
# o Now using getSegments(fit) instead of fit$output.
# 2011-07-10 [HB]
# o Updated code to work with the new column names in PSCBS v0.11.0.
# 2010-10-10 [HB]
# o Added fitC1C2Peaks().
# o Added calibrateC1C2() for PairedPSCBS.
##############################################################################

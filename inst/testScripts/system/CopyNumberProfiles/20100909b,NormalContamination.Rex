library("aroma.cn"); 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAxisBar <- function(x0, x1=NULL, y=par("usr")[3], lwd=4, col=par("col"), length=0.04, xpd=FALSE, ...) {
  if (is.null(x1)) {
    x1 <- x0[2];
    x0 <- x0[1];
  }
  arrows(x0=x0, x1=x1, y0=y, code=3, angle=90, length=length, lwd=lwd, col=col, xpd=xpd, lend=2, ...);
} # xAxisBar()

# Get apply
gapply <- function(X, ..., FUN="[[", useNames=FALSE) {
  res <- sapply(X=X, FUN=FUN, ...);
  if (!useNames) {
    names(res) <- NULL;
  }
  res;
} # gapply()


verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup normal and tumor profiles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The normal/germline CN profile (TCN=2 everywhere)
pN <- cnr(1,1000, 2);
print(pN);

# Build up tumor CN profile generation by generation
pT <- pN + 
      cnr(400,500) + 
      cnr(600,800) +
      cnr(600,700) + 
      cnr(100,200) - cnr(850,900);
print(pT);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnN <- simulateRawCopyNumbers(pN, n=2000, sd=1/3);
cnT <- simulateRawCopyNumbers(pT, n=2000, sd=1/3);


profiles <- list(
  "normal"=list(cnr=pN, cn=cnN, col="black", lty=3),
  "100% tumor"=list(cnr=pT, cn=cnT, col="blue", lty=3)
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build contaminated data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
kappas <- c(0.2, 0.4, 0.8);
cols <- c("magenta", "red", "orange");
ltys <- c(1,1,1);
for (kk in seq(along=kappas)) {
  kappa <- kappas[kk];
  key <- sprintf("%d%% tumor", round(100*(1-kappa)));
  cnr <- kappa * pN + (1-kappa)*pT;
  cn <- kappa * cnN + (1-kappa)*cnT;
  profiles[[key]] <- list(cnr=cnr, cn=cn, col=cols[kk], lty=ltys[kk]);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot profiles with data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subplots(length(profiles), ncol=1);
par(mar=c(0.1,4,2,2)+0.1);
for (kk in seq(along=profiles)) {
  p <- profiles[[kk]];
  xlim <- c(-50, xMax(p$cn)/1e6);
  plot(p$cn, col="gray", xlim=xlim, ylim=c(0,5), axes=FALSE);
  axis(side=2);
  lines(p$cnr, col="white", lty=1, lwd=8);
  lines(p$cnr, col=p$col, lty=1, lwd=5);
  if (kk >= 3) {
    lines(pT, col="blue", lty=3, lwd=3);
  }

  drawDensity(p$cnr, col=p$col, height=50, offset=5);

  title(main=names(profiles)[kk]);
} # for (kk ...)


## # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## # Plot profiles
## # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## subplots(length(profiles)-1, ncol=1);
## par(mar=c(0.1,4,2,2)+0.1);
## for (kk in 2:length(profiles)) {
##   ps <- profiles[unique(c(1,2,kk))];
##   plot(ps[[1]]$cnr, col=ps[[1]]$col, lty=ps[[1]]$lty, lwd=5, axes=FALSE);
##   axis(side=2);
##   lapply(ps[-1], FUN=function(p) lines(p$cnr, col=p$col, lty=p$lty, lwd=5));
##   legend("topright", names(ps), col=gapply(ps,"col"), lty=gapply(ps,"lty"), lwd=5, bty="n");
## } # for (kk ...)


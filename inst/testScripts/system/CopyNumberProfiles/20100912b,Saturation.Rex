library("aroma.cn"); 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
saturate <- function(x, gamma=3) {
  xMax <- max(x, na.rm=TRUE);
  x <- x / xMax;
  y <- gamma * (pnorm(x)-1/2)/(pnorm(gamma)-1/2); 
  y <- y * xMax;
  y;
} # saturate()


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


pS <- clone(pT);
pS$mean <- saturate(pS$mean, gamma=gamma);
pS$mean <- 2 * pS$mean / median(pS$mean, na.rm=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnT <- simulateRawCopyNumbers(pT, n=2000, sd=1/3);
cnS <- simulateRawCopyNumbers(pS, n=2000, sd=1/3);
profiles <- list(
  "Tumor"=list(cnr=pT, cn=cnT, col="blue", lty=3),
  "Tumor w/ saturation"=list(cnr=pS, cn=cnS, col="red", lty=3)
);


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
  lines(pT, col="blue", lty=3, lwd=3);

  drawDensity(p$cnr, col=p$col, height=50, offset=5);

  title(main=names(profiles)[kk]);
} # for (kk ...)
